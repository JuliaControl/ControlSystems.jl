import OrdinaryDiffEq: ODEProblem, Tsit5, DiffEqBase, solve, BS3
import ControlSystemsBase: lsim, step, impulse, HammersteinWienerSystem, DelayLtiSystem, PartitionedStateSpace, SimResult
import DelayDiffEq: MethodOfSteps
# Function for DifferentialEquations lsim
"""
    f_lsim(dx, x, p, t)

Internal function: Dynamics equation for simulation of a linear system.

# Arguments:
- `dx`: State derivative vector written to in place.
- `x`: State
- `p`: is equal to `(A, B, u)` where `u(x, t)` returns the control input
- `t`: Time
"""
@inline function f_lsim(dx, x, p, t) 
    A, B, u = p
    # dx .= A * x .+ B * u(x, t)
    mul!(dx, A, x)
    mul!(dx, B, u(x, t), 1, 1)
end

# This method is more specific than the lsim in ControlSystemsBase that does not specify Continuous timeevol for sys, hence, if ControlSystems is loaded, ControlSystems.lsim will take precedence over ControlSystemsBase.lsim for Continuous systems
function lsim(sys::AbstractStateSpace{Continuous}, u::Function, t::AbstractVector;
        x0::AbstractVecOrMat=zeros(Bool, nstates(sys)), method::Symbol=:cont, alg = Tsit5(), kwargs...)
    ny, nu = size(sys)
    nx = sys.nx
    u0 = u(x0,t[1])
    length(x0) == nx ||
        error("x0 must have length $nx: got length $(length(x0))")
    isa(u0, Number) &&
        error("u must be a vector of size ($nu,)")
    size(u0) == (nu,) ||
        error("u must have size ($nu,): got size $(size(u0))")

    T = promote_type(Float64, eltype(x0), numeric_type(sys))

    dt = t[2] - t[1]

    if method === :zoh
        simsys = c2d(sys, dt, :zoh)
        x,uout = ControlSystemsBase.ltitr(simsys.A, simsys.B, u, t, T.(x0))
    else
        p = (sys.A, sys.B, u)
        sol = solve(ODEProblem(f_lsim, x0, (t[1], t[end]+dt/2), p), alg; saveat=t, kwargs...)
        x = reduce(hcat, sol.u)::Matrix{T}
        uout = Matrix{T}(undef, nu, length(t))
        for i = eachindex(t)
            uout[:, i] = u(@view(x[:, i]), t[i])
        end
        simsys = sys
    end
    y = sys.C*x
    if !iszero(sys.D)
        mul!(y, sys.D, uout, 1, 1)
    end
    return SimResult(y, t, x, uout, simsys) # saves the system that actually produced the simulation
end


## DelayLtiSystem


"""
    res = lsim(sys::DelayLtiSystem, u, t::AbstractArray{<:Real}; x0=fill(0.0, nstates(sys)), alg=MethodOfSteps(Tsit5()), abstol=1e-6, reltol=1e-6, force_dtmin=true, kwargs...)

Simulate system `sys`, over time `t`, using input signal `u`, with initial state `x0`, using method `alg` .

Arguments:

`t`: Has to be an `AbstractVector` with equidistant time samples (`t[i] - t[i-1]` constant)
`u`: Function to determine control signal `ut` at a time `t`, on any of the following forms:
- `u`: Function to determine control signal `uₜ` at a time `t`, on any of the following forms:
    - A constant `Number` or `Vector`, interpreted as a constant input.
    - Function `u(x, t)` that takes the internal state and time, note, the state representation for delay systems is not the same as for rational systems.
    - In-place function `u(uₜ, x, t)`. (Slightly more effienct)

`alg, abstol, reltol` and `kwargs...`: are sent to `DelayDiffEq.solve`.

This methods sets `force_dtmin=true` by default to handle the discontinuity implied by, e.g., step inputs. This may lead to the solver taking a long time to solve ill-conditioned problems rather than exiting early with a warning.

Returns an instance of [`SimResult`](@ref) which can be plotted directly or destructured into `y, t, x, u = res`.
"""
function lsim(sys::DelayLtiSystem{T,S}, u, t::AbstractArray{<:Real};
        x0=fill(zero(T), nstates(sys)),
        alg=DelayDiffEq.MethodOfSteps(Tsit5()),
        abstol=1e-6, reltol=1e-6, force_dtmin = true,
        kwargs...) where {T,S}

    # Make u! in-place function of u
    u! = if isa(u, Number) || isa(u,AbstractVector) # Allow for u to be a constant number or vector
        (isa(u,AbstractVector) && length(u) == sys.nu1) || error("u must have size ($sys.nu1,1): got size $(size(u))")
        let u = u
            @inline (uout, x, t) -> copyto!(uout, u)
        end
    elseif DiffEqBase.isinplace(u, 3) # If u is an inplace (more than 2 argument function)
        u
    else
        let u = u
            @inline (out, x, t) -> copyto!(out, u(x, t))
        end
    end

    _lsim(sys, u!, t, x0, alg; abstol, reltol, force_dtmin, kwargs...)
end

# Generic parametrized dde used for simulating DelayLtiSystem
# We simulate the integral of the states
# u, du are from the notation of variables in DifferentialEquations.jl
# The state from the control system is x
# u!, uout is the control law and its output
function dde_param(du, u, h, p, t)
    A, B1, B2, C1, C2, D11, D12, D21, D22, y, Tau, u!, uout, hout, tmpy, tsave = p

    nx = size(A,1)
    nd = length(Tau)
    ny = size(C1,1)
    
    dx = view(du, 1:nx)
    dY = view(du, (nx+1):(nx+ny))
    dD = view(du, (nx+ny+1):(nx+ny+nd))
    x = view(u, 1:nx)

    # uout = u(t)
    u!(uout, x, t)

    # hout = d(t-Tau)
    for k=1:length(Tau)
        # Get the derivative of the history for d at indices corresponding to Tau
        hout[k] = h(p, t-Tau[k], Val{1}, idxs=(nx+ny+k))
    end

    #dx(t) .= A*x(t) + B1*u(t) +B2*d(t-tau)
    mul!(dx, A, x)
    mul!(dx, B1, uout, true, true)
    mul!(dx, B2, hout, true, true)

    # dY = y(t) = C1*x + D11*u(t) + D12*d(t-Tau)
    mul!(dY, C1, x)
    mul!(dY, D11, uout, true, true)
    mul!(dY, D12, hout, true, true)

    # dD = d(t) = C2*x + D21*u(t) + D22*d(t-Tau)
    mul!(dD, C2, x)
    mul!(dD, D21, uout, true, true)
    mul!(dD, D22, hout, true, true)

    # Save y value in tmpy to be used by dde_saver
    if t in tsave
        # The value of y at this time is given by the derivative
        tmpy .= dY
    end
    return
end

# Save x(t) and y(t) to output
function dde_saver(u,t,integrator)
    A, B1, B2, C1, C2, D11, D12, D21, D22, y, Tau, u!, uout, hout, tmpy, tsave = integrator.p
    nx = size(A,1)
    # y is already saved in tmpy
    u[1:nx], copy(tmpy)
end

function _lsim(sys::DelayLtiSystem{T,S}, Base.@nospecialize(u!), t::AbstractArray{<:Real}, x0::Vector{T}, alg; kwargs...) where {T,S}

    P = sys.P

    t0 = first(t)
    dt = t[2] - t[1]

    # Get all matrices to save on allocations
    A, B1, B2, C1, C2, D11, D12, D21, D22 = P.A, P.B1, P.B2, P.C1, P.C2, P.D11, P.D12, P.D21, P.D22
    Tau = sys.Tau
    nx = size(A,1)
    nd = length(Tau)
    ny = size(C1,1)
    nu = size(B1,2)
    nt = length(t)

    hout = fill(zero(T), nd)            # in place storage for delays
    uout = fill(zero(T), ninputs(sys))  # in place storage for u
    tmpy = similar(x0, ny)              # in place storage for output
    y = Matrix{T}(undef, ny, nt)        # Output matrix
    x = Matrix{T}(undef, nx, nt)        # State matrix

    p = (A, B1, B2, C1, C2, D11, D12, D21, D22, y, Tau, u!, uout, hout, tmpy, t)

    # This callback computes and stores the delay term
    sv = SavedValues(eltype(t), Tuple{Vector{T},Vector{T}})
    cb = SavingCallback(dde_saver, sv, saveat = t)

    # History function, only used for d
    # The true d(t) is stored as a derivative
    h!(p, t, deriv::Type{Val{1}}; idxs=0) = zero(T)

    # The states are x(t), Y(t), D(t), where Y, D are integrals of y(t), d(t)
    u0 = [x0;zeros(T,ny);zeros(T,nd)]
    prob = DelayDiffEq.DDEProblem{true}(dde_param, u0, h!,
                (T(t[1]), T(t[end])),
                p,
                constant_lags=sort(Tau),# Not sure if sort needed
                neutral=true,           # We have derivatives on RHS (d(t)) 
                callback=cb)
    # Important to stop at t since we can not access derivatives in SavingCallback
    sol = DelayDiffEq.solve(prob, alg; tstops=t, saveat=t, kwargs...)

    # Retrive the saved values
    uout2 = zeros(T, nu, nt)
    for k = 1:nt
        x[:,k] .= sv.saveval[k][1]
        y[:,k] .= sv.saveval[k][2]
        @views u!(uout2[:, k], x[:, k], t[k])
    end

    return SimResult(y, t, x, uout2, sys)
end

function Base.step(sys::DelayLtiSystem{T}, t::AbstractVector; kwargs...) where T
    nu = ninputs(sys)
    if t[1] != 0
        throw(ArgumentError("First time point must be 0 in step"))
    end
    u = (out, x, t) -> (t < 0 ? out .= 0 : out .= 1)
    x0=fill(zero(T), nstates(sys))
    if nu == 1
        return lsim(sys, u, t; x0=x0, kwargs...)
    else
        x = Array{T}(undef, nstates(sys), length(t), nu)
        y = Array{T}(undef, noutputs(sys), length(t), nu)
        uout = zeros(T, ninputs(sys), length(t), nu)
        for i=1:nu
            y[:,:,i], tout, x[:,:,i], uout[i,:,i] = lsim(sys[:,i], u, t; x0=x0, kwargs...)
        end
        return SimResult(y, t, x, uout, sys)
    end
end


function impulse(sys::DelayLtiSystem{T}, t::AbstractVector; alg=MethodOfSteps(BS3()), kwargs...) where T
    nu = ninputs(sys)
    iszero(sys.P.D22) || @warn "Impulse with a direct term from delay vector to delay vector can lead to poor results." maxlog=10
    iszero(sys.P.D21) || throw(ArgumentError("Impulse with a direct term from input to delay vector is not implemented. Move the delays to the output instead of input if possible."))
    if t[1] != 0
        throw(ArgumentError("First time point must be 0 in impulse"))
    end
    u = (out, x, t) -> (out .= 0)
    if nu == 1
        return lsim(sys, u, t; alg=alg, x0=sys.P.B[:,1], kwargs...)
    else
        x = Array{T}(undef, nstates(sys), length(t), nu)
        y = Array{T}(undef, noutputs(sys), length(t), nu)
        uout = zeros(T, ninputs(sys), length(t), nu)
        for i=1:nu
            y[:,:,i], tout, x[:,:,i], uout[i,:,i] = lsim(sys[:,i], u, t; alg=alg, x0=sys.P.B[:,i], kwargs...)
        end
        SimResult(y, t, x, uout, sys)
    end
end


## HammersteinWiener


_equation_order_error(D22) = error("D22 contains algebraic loops. This problem can appear if systems with direct feedthrough are used together with nonlinearities on both inputs and outputs. If the offending system is a controller with direct feedthrough, try adding an additional low-pass pole. D22 = $D22")

"""
    _equation_order(D22)

Internal function. Return a permutation vector that determines the order of evaluation such that each equation is evaluated after all of its dependencies.
"""
function _equation_order(D22)
    # NOTE: I think it would theoretically be possible to handle some kinds of algebraic loops, as long as the sub-matrix of the connected component is invertible. Consider the case feedback(ss(1)) == 0.5 which is clearly solvable, while feedback(nonlinearity(identity)) will detect and fail on an algebraic loop. Handling this case may require solving nonlinear systems during simulation and might come with a considerable complexity cost.
    n = LinearAlgebra.checksquare(D22)
    n == 0 && return Int[] # no equations to solve, let's go to the pub and have a mid-strength beer
    all(iszero, D22[diagind(D22)]) || _equation_order_error(D22)

    deps = _get_eq_deps.(Ref(D22), 1:n) # get dependencies of all equations
    # at least one equation must be without dependencies, if not, there is an algebraic loop
    i = findfirst(isempty, deps)
    i === nothing && _equation_order_error(D22)
    order = [i] # all equations present in order can now be considered solved
    for iter = 1:n # should require at most n iterations
        # look for an equation that can be solved given the equations already solved in order
        for j in 1:n
            j ∈ order && continue # already done
            depsj = deps[j]
            if isempty(depsj) || all(d ∈ order for d in depsj)
                # Either no dependencies and can be solved at any time,
                # or we can solve this equation since all of it's dependencies are already solved
                push!(order, j) 
                continue # move on to look at next equation
            end
        end
        if length(order) == n
            @assert sort(order) == 1:n # verify that order is a permutation
            return order # we have solved all equations
        end
    end
    _equation_order_error(D22) # we failed to solve some equations due to algebraic loop
end

_get_eq_deps(D22, i) = findall(!=(0), D22[i,:])

"""
    lsim(sys::HammersteinWienerSystem, u, t::AbstractArray{<:Real}; x0=fill(0.0, nstates(sys)), alg=Tsit5(), abstol=1e-6, reltol=1e-6, kwargs...)

Simulate system `sys`, over time `t`, using input signal `u`, with initial state `x0`, using method `alg` .

# Arguments:

- `t`: Has to be an `AbstractVector` with equidistant time samples (`t[i] - t[i-1]` constant)
- `u`: Function to determine control signal `uₜ` at a time `t`, on any of the following forms:
    Can be a constant `Number` or `Vector`, interpreted as `uₜ .= u` , or
    Function `uₜ .= u(x, t)`, or
    In-place function `u(uₜ, x, t)`. (Slightly more effienct)
- `alg, abstol, reltol` and `kwargs...`: are sent to `OrdinaryDiffEq.solve`.

Returns an instance of [`SimResult`](@ref).
"""
function lsim(sys::HammersteinWienerSystem{T}, u, t::AbstractArray{<:Real};
        x0=fill(zero(T), nstates(sys)),
        alg=Tsit5(),
        abstol=1e-6, reltol=1e-6,
        kwargs...) where {T}

    # Make u! in-place function of u
    u! = if isa(u, Number) || isa(u,AbstractVector) # Allow for u to be a constant number or vector
        (isa(u,AbstractVector) && length(u) == sys.nu1) || error("Vector u must be of length $(sys.nu1)")
        let u = u
            @inline (uout, x, t) -> copyto!(uout, u)
        end
    elseif DiffEqBase.isinplace(u, 3) # If u is an inplace (more than 2 argument function)
        u
    else                              # If u is a regular u(t) function
        let u = u
            @inline (out, x, t) -> copyto!(out, u(x, t))
        end
    end

    _lsim(sys, u!, t, x0, alg; abstol=abstol, reltol=reltol, kwargs...)
end

"Internal function. The right-hand side dynamics function for simulation of HammersteinWienerSystem"
function hw_f(dx, x, p, t)
    (A, B1, B2, C1, C2, D11, D12, D21, D22, f, u!, uout, Δy, order) = p
    u!(uout, x, t)

    # Δy = f = C2*x + D21*u(t) + D22*f(uc)
    mul!(Δy, C2, x)
    mul!(Δy, D21, uout, true, true)
    Δu = Δy # change name, we reuse Δy storage even if the signal after the nonlinearity is Δu
    nonlinear_activation!(Δu, f, D22, order)

    mul!(dx, A, x)
    mul!(dx, B1, uout, true, true)
    mul!(dx, B2, Δu, true, true)
end

function nonlinear_activation!(Δu, f, D22, order)
    n = LinearAlgebra.checksquare(D22)
    length(f) == length(order) == size(Δu, 1) || throw(ArgumentError("inconsistent size of inputs"))
    @inbounds if iszero(D22)
        for k in axes(Δu, 2), i in axes(Δu, 1)
            Δu[i,k] = f[i](Δu[i,k])
        end
    else
        # If D22 has entries on the diagonal, there would be an algebraic loop, 
        # D22[:,order] will however always be lower triangular due to s1*s2 for PartitionedStateSpace, we can thus use the permutation as the order in which to apply nonlinearities
        for k in axes(Δu, 2), i in order 
            Δu[i,k] = f[i](Δu[i,k])
            for j = 1:n
                Δu[j,k] += D22[j, i]*Δu[i,k]
            end
        end
    end
end

function _lsim(sys::HammersteinWienerSystem{T}, u!, t::AbstractArray{<:Real}, x0::Vector{T}, alg; kwargs...) where {T}

    P = sys.P

    t0 = first(t)
    dt = t[2] - t[1]

    # Get all matrices to save on allocations
    A, B1, B2, C1, C2, D11, D12, D21, D22 = P.A, P.B1, P.B2, P.C1, P.C2, P.D11, P.D12, P.D21, P.D22
    f = (sys.f...,) # Turn f into a tuple so that the compiler can specialize on the functions. It's important to not do this before simulation since then there will be a lot of compilation for each unique system.
    nx = size(A,1)
    nd = length(f)
    ny = size(C1,1)
    nu = size(B1,2)
    nt = length(t)

    dy   = fill(zero(T), nd)            # in place storage for nonlinear activations
    uout = fill(zero(T), ninputs(sys))  # in place storage for u
    order = _equation_order(D22) # The order in which to apply nonlinearities

    if nx > 0
        p = (A, B1, B2, C1, C2, D11, D12, D21, D22, f, u!, uout, dy, order)
        prob = ODEProblem{true}(hw_f, x0, (T(t[1]), T(t[end]+dt/2)), p)
        sol = OrdinaryDiffEq.solve(prob, alg; saveat=t, kwargs...)
        x = reduce(hcat, sol.u)::Matrix{T}
    else
        x = zeros(T, 0, nt)        # Empty State matrix
    end
    uout2 = Matrix{T}(undef, nu, length(t))
    for i = eachindex(t)
        @views u!(uout2[:, i], x[:, i], t[i])
    end
    Δy = C2*x
    yout = C1*x
    iszero(D11) || mul!(yout, D11, uout2, true, true)
    
    if !iszero(D12)
        iszero(D21) || mul!(Δy, D21, uout2, true, true)
        nonlinear_activation!(Δy, f, D22, order)
        mul!(yout, D12, Δy, true, true)
    end

    return SimResult(yout, t, x, uout2, sys)
end
