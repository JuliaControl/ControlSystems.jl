# Functions for calculating time response of a system

# XXX : `step` is a function in Base, with a different meaning than it has
# here. This shouldn't be an issue, but it might be.
"""
    y, t, x = step(sys[, tfinal])
    y, t, x = step(sys[, t])

Calculate the step response of system `sys`. If the final time `tfinal` or time
vector `t` is not provided, one is calculated based on the system pole
locations.

`y` has size `(length(t), ny, nu)`, `x` has size `(length(t), nx, nu)`"""
function Base.step(sys::AbstractStateSpace, t::AbstractVector; method=:cont)
    lt = length(t)
    ny, nu = size(sys)
    nx = sys.nx
    u = (x,t)->[one(eltype(t))]
    x0 = zeros(nx)
    if nu == 1
        y, tout, x, _ = lsim(sys, u, t, x0=x0, method=method)
    else
        x = Array{Float64}(undef, lt, nx, nu)
        y = Array{Float64}(undef, lt, ny, nu)
        for i=1:nu
            y[:,:,i], tout, x[:,:,i],_ = lsim(sys[:,i], u, t, x0=x0, method=method)
        end
    end
    return y, t, x
end

Base.step(sys::LTISystem, tfinal::Real; kwargs...) = step(sys, _default_time_vector(sys, tfinal); kwargs...)
Base.step(sys::LTISystem; kwargs...) = step(sys, _default_time_vector(sys); kwargs...)
Base.step(sys::TransferFunction, t::AbstractVector; kwargs...) = step(ss(sys), t::AbstractVector; kwargs...)

"""
    y, t, x = impulse(sys[, tfinal])
    y, t, x = impulse(sys[, t])

Calculate the impulse response of system `sys`. If the final time `tfinal` or time
vector `t` is not provided, one is calculated based on the system pole
locations.

`y` has size `(length(t), ny, nu)`, `x` has size `(length(t), nx, nu)`"""
function impulse(sys::AbstractStateSpace, t::AbstractVector; method=:cont)
    T = promote_type(eltype(sys.A), Float64)
    lt = length(t)
    ny, nu = size(sys)
    nx = sys.nx
    if iscontinuous(sys) #&& method === :cont
        u = (x,t) -> [zero(T)]
        # impulse response equivalent to unforced response of
        # ss(A, 0, C, 0) with x0 = B.
        imp_sys = ss(sys.A, zeros(T, nx, 1), sys.C, zeros(T, ny, 1))
        x0s = sys.B
    else
        u = (x,i) -> (i == t[1] ? [one(T)]/sys.Ts : [zero(T)])
        imp_sys = sys
        x0s = zeros(T, nx, nu)
    end
    if nu == 1 # Why two cases # QUESTION: Not type stable?
        y, t, x,_ = lsim(sys, u, t, x0=x0s[:], method=method)
    else
        x = Array{T}(undef, lt, nx, nu)
        y = Array{T}(undef, lt, ny, nu)
        for i=1:nu
            y[:,:,i], t, x[:,:,i],_ = lsim(sys[:,i], u, t, x0=x0s[:,i], method=method)
        end
    end
    return y, t, x
end

impulse(sys::LTISystem, tfinal::Real; kwags...) = impulse(sys, _default_time_vector(sys, tfinal); kwags...)
impulse(sys::LTISystem; kwags...) = impulse(sys, _default_time_vector(sys); kwags...)
impulse(sys::TransferFunction, t::AbstractVector; kwags...) = impulse(ss(sys), t; kwags...)

"""
    y, t, x = lsim(sys, u[, t]; x0, method])
    y, t, x, uout = lsim(sys, u::Function, t; x0, method)

Calculate the time response of system `sys` to input `u`. If `x0` is ommitted,
a zero vector is used.

`y`, `x`, `uout` has time in the first dimension. Initial state `x0` defaults to zero.

Continuous time systems are simulated using an ODE solver if `u` is a function. If `u` is an array, the system is discretized before simulation. For a lower level inteface, see `?Simulator` and `?solve`

`u` can be a function or a matrix/vector of precalculated control signals.
If `u` is a function, then `u(x,i)` (`u(x,t)`) is called to calculate the control signal every iteration (time instance used by solver). This can be used to provide a control law such as state feedback `u(x,t) = -L*x` calculated by `lqr`.
To simulate a unit step, use `(x,i)-> 1`, for a ramp, use `(x,i)-> i*Ts`, for a step at `t=5`, use (x,i)-> (i*Ts >= 5) etc.

Usage example:
```julia
using LinearAlgebra # For identity matrix I
using Plots

A = [0 1; 0 0]
B = [0;1]
C = [1 0]
sys = ss(A,B,C,0)
Q = I
R = I
L = lqr(sys,Q,R)

u(x,t) = -L*x # Form control law,
t=0:0.1:5
x0 = [1,0]
y, t, x, uout = lsim(sys,u,t,x0=x0)
plot(t,x, lab=["Position" "Velocity"], xlabel="Time [s]")
```
"""
function lsim(sys::AbstractStateSpace, u::AbstractVecOrMat, t::AbstractVector;
        x0::AbstractVector=zeros(Bool, sys.nx), method::Symbol=:unspecified)
    ny, nu = size(sys)
    nx = sys.nx

    if length(x0) != nx
        error("size(x0) must match the number of states of sys")
    end
    if !(size(u) in [(length(t), nu) (length(t),)])
        error("u must be of size (length(t), nu)")
    end

    dt = Float64(t[2] - t[1])
    if !all(x -> x ≈ dt, diff(t))
        error("time vector t must be uniformly spaced")
    end

    if iscontinuous(sys)
        if method === :unspecified
            method = _issmooth(u) ? :foh : :zoh
        end

        if method === :zoh
            dsys = c2d(sys, dt, :zoh)
        elseif method === :foh
            dsys, x0map = c2d_x0map(sys, dt, :foh)
            x0 = x0map*[x0; transpose(u)[:,1]]
        else
            error("Unsupported discretization method")
        end
    else
        if sys.Ts != dt
            error("Time vector must match sample time of discrete-time system")
        end
        dsys = sys
    end

    x = ltitr(dsys.A, dsys.B, u, x0)
    y = transpose(sys.C*transpose(x) + sys.D*transpose(u))
    return y, t, x
end

function lsim(sys::AbstractStateSpace{<:Discrete}, u::AbstractVecOrMat; kwargs...)
    t = range(0, length=size(u, 1), step=sys.Ts)
    lsim(sys, u, t; kwargs...)
end

@deprecate lsim(sys, u, t, x0) lsim(sys, u, t; x0=x0)
@deprecate lsim(sys, u, t, x0, method) lsim(sys, u, t; x0=x0, method=method)

function lsim(sys::AbstractStateSpace, u::Function, tfinal::Real, args...; kwargs...)
    t = _default_time_vector(sys, tfinal)
    lsim(sys, u, t, args...; kwargs...)
end

function lsim(sys::AbstractStateSpace, u::Function, t::AbstractVector;
        x0::VecOrMat=zeros(sys.nx), method::Symbol=:cont)
    ny, nu = size(sys)
    nx = sys.nx
    u0 = u(x0,1)
    if length(x0) != nx
        error("size(x0) must match the number of states of sys")
    elseif !(u0 isa Number && nu == 1) && (size(u0) != (nu,) && size(u0) != (nu,1))
        error("return value of u must be of size nu")
    end
    T = promote_type(Float64, eltype(x0), numeric_type(sys))

    dt = T(t[2] - t[1])

    if !iscontinuous(sys) || method === :zoh
        if iscontinuous(sys)
            dsys = c2d(sys, dt, :zoh)
        else
            if sys.Ts != dt
                error("Time vector must match sample time for discrete system")
            end
            dsys = sys
        end
        x,uout = ltitr(dsys.A, dsys.B, u, t, T.(x0))
    else
        s = Simulator(sys, u)
        sol = solve(s, T.(x0), (t[1],t[end]), Tsit5())
        x = sol(t)'
        uout = Array{eltype(x)}(undef, length(t), ninputs(sys))
        for (i,ti) in enumerate(t)
            uout[i,:] = u(x[i,:],ti)'
        end
    end
    y = transpose(sys.C*transpose(x) + sys.D*transpose(uout))
    return y, t, x, uout
end


lsim(sys::TransferFunction, u, t, args...; kwargs...) = lsim(ss(sys), u, t, args...; kwargs...)


"""
    ltitr(A, B, u[,x0])
    ltitr(A, B, u::Function, iters[,x0])

Simulate the discrete time system `x[k + 1] = A x[k] + B u[k]`, returning `x`.
If `x0` is not provided, a zero-vector is used.

The type of `x0` determines the matrix structure of the returned result,
e.g, `x0` should prefereably not be a sparse vector.

If `u` is a function, then `u(x,i)` is called to calculate the control signal every iteration. This can be used to provide a control law such as state feedback `u=-Lx` calculated by `lqr`. In this case, an integrer `iters` must be provided that indicates the number of iterations.
"""
@views function ltitr(A::AbstractMatrix, B::AbstractMatrix, u::AbstractVecOrMat,
        x0::AbstractVector=zeros(eltype(A), size(A, 1)))

    T = promote_type(LinearAlgebra.promote_op(LinearAlgebra.matprod, eltype(A), eltype(x0)),
                      LinearAlgebra.promote_op(LinearAlgebra.matprod, eltype(B), eltype(u)))

    n = size(u, 1)

    # Transposing u allows column-wise operations, which apparently is faster.
    ut = transpose(u)

    # Using similar instead of Matrix{T} to allow for CuArrays to be used.
    # This approach is problematic if x0 is sparse for example, but was considered
    # to be good enough for now
    x = similar(x0, T, (length(x0), n))

    x[:,1] .= x0
    mul!(x[:, 2:end], B, transpose(u[1:end-1, :])) # Do all multiplications B*u[:,k] to save view allocations

    tmp = similar(x0, T) # temporary vector for storing A*x[:,k]
    for k=1:n-1
        mul!(tmp, A, x[:,k])
        x[:,k+1] .+= tmp
    end
    return transpose(x)
end

function ltitr(A::AbstractMatrix{T}, B::AbstractMatrix{T}, u::Function, t,
    x0::VecOrMat=zeros(T, size(A, 1))) where T
    iters = length(t)
    x = similar(A, size(A, 1), iters)
    uout = similar(A, size(B, 2), iters)

    for i=1:iters
        x[:,i] = x0
        uout[:,i] .= u(x0,t[i])
        x0 = A * x0 + B * uout[:,i]
    end
    return transpose(x), transpose(uout)
end

# HELPERS:

# TODO: This is a poor heuristic to estimate a "good" time vector to use for
# simulation, in cases when one isn't provided.
function _default_time_vector(sys::LTISystem, tfinal::Real=-1)
    dt = _default_dt(sys)
    if tfinal == -1
        tfinal = 200dt
    end
    return 0:dt:tfinal
end

function _default_dt(sys::LTISystem)
    if isdiscrete(sys)
        return sys.Ts
    elseif all(iszero, pole(sys)) # Static or pure integrators
        return 0.05
    else
        ω0_max = maximum(abs.(pole(sys)))
        dt = round(1/(12*ω0_max), sigdigits=2)
        return dt
    end
end



#TODO a reasonable check
_issmooth(u::Function) = false

# Determine if a signal is "smooth"
function _issmooth(u, thresh::AbstractFloat=0.75)
    u = [zeros(1, size(u, 2)); u]       # Start from 0 signal always
    dist = maximum(u) - minimum(u)
    du = abs.(diff(u, dims=1))
    return !isempty(du) && all(maximum(du) <= thresh*dist)
end
