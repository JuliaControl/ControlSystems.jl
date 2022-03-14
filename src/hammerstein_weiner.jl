function freqresp!(R::Array{T,3}, sys::HammersteinWienerSystem, ω::AbstractVector{W}) where {T, W <: Real}
    all(f isa Offset for f in sys.Tau) && return freqresp!(R, lft(sys.P.P, ss(I(length(sys.Tau)), timeevol(sys))), ω)
    throw(ArgumentError("Frequency response is not defined for HammersteinWienerSystem with nonlinearities. Call linearize to obtain a linearized system"))
end

function evalfr(sys::HammersteinWienerSystem, s)
    all(f isa Offset for f in sys.Tau) && return evalfr(lft(sys.P.P, ss(I(length(sys.Tau)), timeevol(sys))), s)
    throw(ArgumentError("Frequency response is not defined for HammersteinWienerSystem with nonlinearities. Call linearize to obtain a linearized system"))
end

function linearize(sys::HammersteinWienerSystem, Δy)
    f = [ForwardDiff.derivative(f, dy) for (f,dy) in zip(sys.Tau, Δy)]
    lft(sys.P.P, ss(diagm(f), timeevol(sys)))
end


"""
    `y, t, x = lsim(sys::HammersteinWienerSystem, u, t::AbstractArray{<:Real}; x0=fill(0.0, nstates(sys)), alg=MethodOfSteps(Tsit5()), abstol=1e-6, reltol=1e-6, kwargs...)`

    Simulate system `sys`, over time `t`, using input signal `u`, with initial state `x0`, using method `alg` .

    Arguments:

    `t`: Has to be an `AbstractVector` with equidistant time samples (`t[i] - t[i-1]` constant)
    `u`: Function to determine control signal `ut` at a time `t`, on any of the following forms:
        Can be a constant `Number` or `Vector`, interpreted as `ut .= u` , or
        Function `ut .= u(t)`, or
        In-place function `u(ut, t)`. (Slightly more effienct)
    `alg, abstol, reltol` and `kwargs...`: are sent to `DelayDiffEq.solve`.

    Returns: times `t`, and `y` and `x` at those times.
"""
function lsim(sys::HammersteinWienerSystem{T}, u, t::AbstractArray{<:Real};
        x0=fill(zero(T), nstates(sys)),
        alg=Tsit5(),
        abstol=1e-6, reltol=1e-6,
        kwargs...) where {T,S}

    # Make u! in-place function of u
    u! = if isa(u, Number) || isa(u,AbstractVector) # Allow for u to be a constant number or vector
        (uout, t) -> uout .= u
    elseif DiffEqBase.isinplace(u, 2)               # If u is an inplace (more than 1 argument function)
        u
    else                                            # If u is a regular u(t) function
        (out, t) -> (out .= u(t))
    end

    _lsim(sys, u!, t, x0, alg; abstol=abstol, reltol=reltol, kwargs...)
end


function hw_f(du, u, p, t)
    (A, B1, B2, C1, C2, D11, D12, D21, D22, Tau, u!, uout, Δy) = p

    nx = size(A,1)
    nd = length(Tau)
    ny = size(C1,1)

    yinds = (nx+1):(nx+ny)
    dx = view(du, 1:nx)

    x = view(u, 1:nx)
    u!(uout, t)

    # Δy = f = C2*x + D21*u(t) + D22*f(uc)
    mul!(Δy, C2, x)
    mul!(Δy, D21, uout, true, true)
    # If D22 has entries on the diagonal, there would be an algebraic loop, 
    # D22 will however always be upper triangular due to s1*s2 for PartitionedStateSpace, we can thus start from the bottom and apply nonlinearities
    Δu = Δy # change name, we reuse Δy storage even if the signal after the nonlinearity is Δu
    for i in reverse(eachindex(Δy))
        Δu[i] = Tau[i](Δy[i]) #
        for j = i-1:-1:1 
            Δu[j] += D22[j, i]*Δu[i]
        end
    end

    mul!(dx, A, x)
    mul!(dx, B1, uout, true, true)
    mul!(dx, B2, Δu, true, true)
end


function _lsim(sys::HammersteinWienerSystem{T}, Base.@nospecialize(u!), t::AbstractArray{<:Real}, x0::Vector{T}, alg; kwargs...) where {T}

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

    dy   = fill(zero(T), nd)            # in place storage for nonlinear activations
    uout = fill(zero(T), ninputs(sys))  # in place storage for u

    if nx > 0
        p = (A, B1, B2, C1, C2, D11, D12, D21, D22, Tau, u!, uout, dy)
        u0 = x0
        prob = ODEProblem{true}(hw_f, u0, (T(t[1]), T(t[end])), p)
        sol = OrdinaryDiffEq.solve(prob, alg; saveat=t, kwargs...)
        x = reduce(hcat, sol.u)::Matrix{T}
    else
        x = Matrix{T}(undef, nx, nt)        # State matrix
    end
    uout2 = Matrix{T}(undef, nu, length(t))
    for i = eachindex(t)
        @views u!(uout2[:, i], t[i])
    end
    Δy = C2*x
    yout = C1*x
    iszero(D11) || mul!(yout, D11, uout2, true, true)
    
    if !iszero(D12)
        iszero(D21) || mul!(Δy, D21, uout2, true, true)
        if iszero(D22)
            for k in axes(Δy, 2), i in reverse(axes(Δy, 1))
                Δy[i,k] = Tau[i](Δy[i,k])
            end
        else
            iszero(LowerTriangular(D22)) ||
                error("D22 must be strict upper triangular to prevent algebraic loops. It's recommended to build nonlinear systems using the constructor `nonlinearity` to enforce this structure.")
            for k in axes(Δy, 2), i in reverse(axes(Δy, 1))
                Δy[i,k] = Tau[i](Δy[i,k])
                for j = i-1:-1:1 
                    Δy[j,k] += D22[j, i]*Δy[i,k]
                end
            end
        end     
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
    if t[1] != 0
        throw(ArgumentError("First time point must be 0 in step"))
    end
    u = (out, t) -> (t < 0 ? out .= 0 : out .= 1)
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


function impulse(sys::HammersteinWienerSystem{T}, t::AbstractVector; kwargs...) where T
    nu = ninputs(sys)
    iszero(sys.P.D22) || @warn "Impulse with a direct term from delay vector to delay vector can lead to poor results." maxlog=10
    iszero(sys.P.D21) || throw(ArgumentError("Impulse with a direct term from input to delay vector is not implemented. Move the delays to the output instead of input if possible."))
    if t[1] != 0
        throw(ArgumentError("First time point must be 0 in impulse"))
    end
    u = (out, t) -> (out .= 0)
    if nu == 1
        return lsim(sys, u, t; alg=alg, x0=sys.P.B[:,1], kwargs...)
    else
        x = Array{T}(undef, nstates(sys), length(t), nu)
        y = Array{T}(undef, noutputs(sys), length(t), nu)
        uout = zeros(T, ninputs(sys), length(t), nu)
        for i=1:nu
            y[:,:,i], tout, x[:,:,i], uout[:,:,i] = lsim(sys[:,i], u, t; alg=alg, x0=sys.P.B[:,i], kwargs...)
        end
        SimResult(y, t, x, uout, sys)
    end
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

Create a saturating nonlinearity.
"""
saturation(args...) = nonlinearity(Saturation(args...))

struct Offset{T} <: Function
    o::T
end

(s::Offset)(x) = x .+ s.o

"""
    offset(val)

Create a constant-offset nonlinearity `x -> x + val`.
"""
offset(val::Number) = nonlinearity(Offset(val))
offset(v::AbstractVector) = nonlinearity(Offset.(v))
