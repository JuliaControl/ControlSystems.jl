# Functions for calculating time response of a system

# XXX : `step` is a function in Base, with a different meaning than it has
# here. This shouldn't be an issue, but it might be.
"""`y, t, x = step(sys[, Tf])` or `y, t, x = step(sys[, t])`

Calculate the step response of system `sys`. If the final time `Tf` or time
vector `t` is not provided, one is calculated based on the system pole
locations.

`y` has size `(length(t), ny, nu)`, `x` has size `(length(t), nx, nu)`"""
function Base.step(sys::StateSpace, t::AbstractVector; method=:cont)
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

Base.step(sys::LTISystem, Tf::Real; kwargs...) = step(sys, _default_time_vector(sys, Tf); kwargs...)
Base.step(sys::LTISystem; kwargs...) = step(sys, _default_time_vector(sys); kwargs...)
Base.step(sys::TransferFunction, t::AbstractVector; kwargs...) = step(ss(sys), t::AbstractVector; kwargs...)

"""`y, t, x = impulse(sys[, Tf])` or `y, t, x = impulse(sys[, t])`

Calculate the impulse response of system `sys`. If the final time `Tf` or time
vector `t` is not provided, one is calculated based on the system pole
locations.

`y` has size `(length(t), ny, nu)`, `x` has size `(length(t), nx, nu)`"""
function impulse(sys::StateSpace, t::AbstractVector; method=:cont)
    T = promote_type(eltype(sys.A), Float64)
    lt = length(t)
    ny, nu = size(sys)
    nx = sys.nx
    if iscontinuous(sys) #&& method == :cont
        u = (x,i) -> [zero(T)]
        # impulse response equivalent to unforced response of
        # ss(A, 0, C, 0) with x0 = B.
        imp_sys = ss(sys.A, zeros(T, nx, 1), sys.C, zeros(T, ny, 1))
        x0s = sys.B
    else
        u = (x,i) -> i == t[1] ? [one(T)]/sys.Ts : [zero(T)]
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

impulse(sys::LTISystem, Tf::Real; kwags...) = impulse(sys, _default_time_vector(sys, Tf); kwags...)
impulse(sys::LTISystem; kwags...) = impulse(sys, _default_time_vector(sys); kwags...)
impulse(sys::TransferFunction, t::AbstractVector; kwags...) = impulse(ss(sys), t; kwags...)

"""`y, t, x = lsim(sys, u, t; x0, method])`

`y, t, x, uout = lsim(sys, u::Function, t; x0, method)`

Calculate the time response of system `sys` to input `u`. If `x0` is ommitted,
a zero vector is used.

`y`, `x`, `uout` has time in the first dimension. Initial state `x0` defaults to zero.

Continuous time systems are simulated using an ODE solver if `u` is a function. If `u` is an array, the system is discretized before simulation. For a lower level inteface, see `?Simulator` and `?solve`

`u` can be a function or a matrix/vector of precalculated control signals.
If `u` is a function, then `u(x,i)` (`u(x,t)`) is called to calculate the control signal every iteration (time instance used by solver). This can be used to provide a control law such as state feedback `u(x,t) = -L*x` calculated by `lqr`.
To simulate a unit step, use `(x,i)-> 1`, for a ramp, use `(x,i)-> i*h`, for a step at `t=5`, use (x,i)-> (i*h >= 5) etc.

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
function lsim(sys::StateSpace, u::AbstractVecOrMat, t::AbstractVector;
        x0::VecOrMat=zeros(sys.nx), method::Symbol=_issmooth(u) ? :foh : :zoh)
    ny, nu = size(sys)
    nx = sys.nx

    if length(x0) != nx
        error("size(x0) must match the number of states of sys")
    end
    if !(size(u) in [(length(t), nu) (length(t),)])
        error("u must be of size (length(t), nu)")
    end

    dt = Float64(t[2] - t[1])
    if !iscontinuous(sys) || method == :zoh
        if iscontinuous(sys) # Looks strange to check iscontinuous again
            dsys = c2d(sys, dt, :zoh)[1]
        else
            if sys.Ts != dt
                error("Time vector must match sample time for discrete system")
            end
            dsys = sys
        end
    else
        dsys, x0map = c2d(sys, dt, :foh)
        x0 = x0map*[x0; transpose(u[1:1,:])]
    end
    x = ltitr(dsys.A, dsys.B, Float64.(u), Float64.(x0))
    y = transpose(sys.C*transpose(x) + sys.D*transpose(u))
    return y, t, x
end

@deprecate lsim(sys, u, t, x0) lsim(sys, u, t; x0=x0)
@deprecate lsim(sys, u, t, x0, method) lsim(sys, u, t; x0=x0, method=method)

function lsim(sys::StateSpace, u::Function, t::AbstractVector;
        x0::VecOrMat=zeros(sys.nx), method::Symbol=:cont)
    ny, nu = size(sys)
    nx = sys.nx
    if length(x0) != nx
        error("size(x0) must match the number of states of sys")
    elseif size(u(x0,1)) != (nu,) && size(u(x0,1)) != (nu,1)
        error("return value of u must be of size nu")
    end
    T = promote_type(Float64, eltype(x0))

    dt = T(t[2] - t[1])
    if !iscontinuous(sys) || method == :zoh
        if iscontinuous(sys)
            dsys = c2d(sys, dt, :zoh)[1]
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


"""`ltitr(A, B, u[,x0])`

`ltitr(A, B, u::Function, iters[,x0])`

Simulate the discrete time system `x[k + 1] = A x[k] + B u[k]`, returning `x`.
If `x0` is not provided, a zero-vector is used.

If `u` is a function, then `u(x,i)` is called to calculate the control signal every iteration. This can be used to provide a control law such as state feedback `u=-Lx` calculated by `lqr`. In this case, an integrer `iters` must be provided that indicates the number of iterations.
"""
function ltitr(A::AbstractMatrix{T}, B::AbstractMatrix{T}, u::AbstractVecOrMat,
        x0::VecOrMat=zeros(T, size(A, 1))) where T
    n = size(u, 1)
    x = similar(A, size(A, 1), n)
    for i=1:n
        x[:,i] = x0
        x0 = A * x0 + B * u[i,:]
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
        uout[:,i] = u(x0,t[i])
        x0 = A * x0 + B * uout[:,i]
    end
    return transpose(x), transpose(uout)
end

# HELPERS:

# TODO: This is a poor heuristic to estimate a "good" time vector to use for
# simulation, in cases when one isn't provided.
function _default_time_vector(sys::LTISystem, Tf::Real=-1)
    Ts = _default_Ts(sys)
    if Tf == -1
        Tf = 100*Ts
    end
    return 0:Ts:Tf
end

function _default_Ts(sys::LTISystem)
    if !iscontinuous(sys)
        Ts = sys.Ts
    elseif !isstable(sys)
        Ts = 0.05
    else
        ps = pole(sys)
        r = minimum([abs.(real.(ps));0])
        if r == 0.0
            r = 1.0
        end
        Ts = 0.07/r
    end
    return Ts
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
