function freqresp!(R::Array{T,3}, sys::HammersteinWienerSystem, ω::AbstractVector{W}) where {T, W <: Real}
    all(f isa Offset for f in sys.f) && return freqresp!(R, lft(sys.P.P, ss(I(length(sys.f)), timeevol(sys))), ω)
    throw(ArgumentError("Frequency response is not defined for HammersteinWienerSystem with nonlinearities. Call linearize to obtain a linearized system"))
end

function evalfr(sys::HammersteinWienerSystem, s)
    all(f isa Offset for f in sys.f) && return evalfr(lft(sys.P.P, ss(I(length(sys.f)), timeevol(sys))), s)
    throw(ArgumentError("Frequency response is not defined for HammersteinWienerSystem with nonlinearities. Call linearize to obtain a linearized system"))
end

"""
    linearize(sys::HammersteinWienerSystem, Δy)

Linearize the nonlinear system `sys` around the operating point implied by the specified Δy
```
      ┌─────────┐
 y◄───┤         │◄────u
      │    P    │
Δy┌───┤         │◄───┐Δu
  │   └─────────┘    │
  │                  │
  │      ┌───┐       │
  │      │   │       │
  └─────►│ f ├───────┘
         │   │
         └───┘
```
$nonlinear_warning
"""
function linearize(sys::HammersteinWienerSystem, Δy)
    f = [ForwardDiff.derivative(f, dy) for (f,dy) in zip(sys.f, Δy)]
    lft(sys.P.P, ss(diagm(f), timeevol(sys)))
end

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

function _bounds_and_features(sys::HammersteinWienerSystem, plot::Val)
    _bounds_and_features(sys.P.P, plot)
end

# Againm we have to do something for default vectors, more or less a copy from timeresp.jl
function _default_dt(sys::HammersteinWienerSystem)
    _default_dt(sys.P.P)
end

function Base.step(sys::HammersteinWienerSystem{T}, t::AbstractVector; kwargs...) where T
    nu = ninputs(sys)
    t[1] == 0 || throw(ArgumentError("First time point must be 0 in step"))
    u = (out, x, t) -> (t < 0 ? out .= 0 : out .= 1)
    x0 = zeros(T, nstates(sys))
    x = Array{T}(undef, nstates(sys), length(t), nu)
    y = Array{T}(undef, noutputs(sys), length(t), nu)
    uout = zeros(T, ninputs(sys), length(t), nu)
    for i=1:nu
        y[:,:,i], tout, x[:,:,i], uout[i,:,i] = lsim(sys[:,i], u, t; x0=x0, kwargs...)
    end
    return SimResult(y, t, x, uout, sys)
end


function impulse(sys::HammersteinWienerSystem{T}, t::AbstractVector; kwargs...) where T
    nu = ninputs(sys)
    t[1] == 0 || throw(ArgumentError("First time point must be 0 in impulse"))
    u = @inline (out, x, t) -> (out .= 0)
    x = Array{T}(undef, nstates(sys), length(t), nu)
    y = Array{T}(undef, noutputs(sys), length(t), nu)
    uout = zeros(T, ninputs(sys), length(t), nu)
    for i=1:nu
        y[:,:,i], tout, x[:,:,i], uout[i,:,i] = lsim(sys[:,i], u, t; x0=sys.P.B[:,i], kwargs...)
    end
    SimResult(y, t, x, uout, sys)
end


struct Saturation{T} <: Function
    l::T
    u::T
end
Saturation(u) = Saturation(-u, u)

(s::Saturation)(x) = clamp(x, s.l, s.u)

"""
    saturation(val)
    saturation(lower, upper)

Create a saturating nonlinearity. Connect it to the output of a controller `C` using
```
Csat = saturation(val) * C
```

$nonlinear_warning
"""
saturation(args...) = nonlinearity(Saturation(args...))
saturation(v::AbstractVector, args...) = nonlinearity(Saturation.(v, args...))
Base.show(io::IO, f::Saturation) = f.u == -f.l ? print(io, "saturation($(f.u))") : print(io, "saturation($(f.l), $(f.u))")

struct Offset{T} <: Function
    o::T
end

(s::Offset)(x) = x .+ s.o

"""
    offset(val)

Create a constant-offset nonlinearity `x -> x + val`.

$nonlinear_warning
"""
offset(val::Number) = nonlinearity(Offset(val))
offset(v::AbstractVector) = nonlinearity(Offset.(v))

Base.show(io::IO, f::Offset) = print(io, "offset($(f.o))")