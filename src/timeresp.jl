# Functions for calculating time response of a system

# XXX : `step` is a function in Base, with a different meaning than it has
# here. This shouldn't be an issue, but it might be.
@doc """`[y, t, x] = step(sys[, Tf])` or `[y, t, x] = step(sys[, t])`

Calculate the step response of system `sys`. If the final time `Tf` or time
vector `t` is not provided, one is calculated based on the system pole
locations.""" ->
function Base.step(sys::StateSpace, t::AbstractVector)
    lt = length(t)
    ny, nu = size(sys)
    nx = sys.nx
    u = ones(size(t))
    x0 = zeros(nx, 1)
    if nu == 1
        y, t, x = lsim(sys, u, t, x0, :zoh)
    else
        x = Array(Float64, lt, nx, nu)
        y = Array(Float64, lt, ny, nu)
        for i=1:nu
            y[:,:,i], t, x[:,:,i] = lsim(sys[:,i], u, t, x0, :zoh)
        end
    end
    return y, t, x
end
function Base.step(sys::TransferFunction{SisoGeneralized}, t::AbstractVector)
    lsim(sys, ones(length(t), sys.nu), t)
end

Base.step(sys::LTISystem, Tf::Real) = step(sys, _default_time_vector(sys, Tf))
Base.step(sys::LTISystem) = step(sys, _default_time_vector(sys))
Base.step(sys::TransferFunction, t::AbstractVector) = step(ss(sys), t::AbstractVector)

@doc """`[y, t, x] = impulse(sys[, Tf])` or `[y, t, x] = impulse(sys[, t])`

Calculate the impulse response of system `sys`. If the final time `Tf` or time
vector `t` is not provided, one is calculated based on the system pole
locations.""" ->
function impulse(sys::StateSpace, t::AbstractVector)
    lt = length(t)
    ny, nu = size(sys)
    nx = sys.nx
    u = zeros(size(t))
    if iscontinuous(sys)
        # impulse response equivalent to unforced response of
        # ss(A, 0, C, 0) with x0 = B.
        imp_sys = ss(sys.A, zeros(nx, 1), sys.C, zeros(ny, 1))
        x0s = sys.B
    else
        imp_sys = sys
        x0s = zeros(nx, nu)
        u[1] = 1/sys.Ts
    end
    if nu == 1
        y, t, x = lsim(sys, u, t, x0s, :zoh)
    else
        x = Array(Float64, lt, nx, nu)
        y = Array(Float64, lt, ny, nu)
        for i=1:nu
            y[:,:,i], t, x[:,:,i] = lsim(sys[:,i], u, t, x0s[:,i], :zoh)
        end
    end
    return y, t, x
end
function impulse(sys::TransferFunction{SisoGeneralized}, t::AbstractVector)
    u = zeros(length(t), sys.nu)
    u[1,:] = 1/(t[2]-t[1])
    lsim(sys::TransferFunction{SisoGeneralized}, u, t)
end

impulse(sys::LTISystem, Tf::Real) = impulse(sys, _default_time_vector(sys, Tf))
impulse(sys::LTISystem) = impulse(sys, _default_time_vector(sys))
impulse(sys::TransferFunction, t::AbstractVector) = impulse(ss(sys), t)

@doc """`y, t, x = lsim(sys, u, t[, x0, method])`

`y, t, x, uout = lsim(sys, u::Function, t[, x0, method])`

Calculate the time response of system `sys` to input `u`. If `x0` is ommitted,
a zero vector is used.

Continuous time systems are discretized before simulation. By default, the
method is chosen based on the smoothness of the input signal. Optionally, the
`method` parameter can be specified as either `:zoh` or `:foh`.

`u` can be a function or a matrix/vector of precalculated control signals.
If `u` is a function, then `u(i,x)` is called to calculate the control signal every iteration. This can be used to provide a control law such as state feedback `u(t,x) = -L*x` calculated by `lqr`.
To simulate a unit step, use `(i,x)-> 1`, for a ramp, use `(i,x)-> i*h`, for a step at `t=5`, use (i,x)-> (i*h >= 5) etc.

Usage example:
```julia
A = [0 1; 0 0]
B = [0;1]
C = [1 0]
sys = ss(A,B,C,0)
Q = eye(2)
R = eye(1)
L = lqr(sys,Q,R)

u(t,x) = -L*x # Form control law,
t=0:0.1:5
x0 = [1,0]
y, t, x, uout = lsim(sys,u,t,x0)
plot(t,x, lab=["Position", "Velocity"]', xlabel="Time [s]")
```
""" ->
function lsim(sys::StateSpace, u::AbstractVecOrMat, t::AbstractVector,
        x0::VecOrMat=zeros(sys.nx, 1), method::Symbol=_issmooth(u) ? :foh : :zoh)
    ny, nu = size(sys)
    nx = sys.nx

    if length(x0) != nx
        error("size(x0) must match the number of states of sys")
    elseif !any(size(u) .== [(length(t), nu) (length(t),)])
        error("u must be of size (length(t), nu)")
    end

    dt = Float64(t[2] - t[1])
    if !iscontinuous(sys) || method == :zoh
        if iscontinuous(sys)
            dsys = c2d(sys, dt, :zoh)[1]
        else
            if sys.Ts != dt
                error("Time vector must match sample time for discrete system")
            end
            dsys = sys
        end
    else
        dsys, x0map = c2d(sys, dt, :foh)
        x0 = x0map*[x0; u[1,:].']
    end
    x = ltitr(dsys.A, dsys.B, map(Float64,u), map(Float64,x0))
    y = (sys.C*(x.') + sys.D*(u.')).'
    return y, t, x
end

function lsim(sys::StateSpace, u::Function, t::AbstractVector,
        x0::VecOrMat=zeros(sys.nx, 1), method::Symbol=:zoh)
    ny, nu = size(sys)
    nx = sys.nx
    if length(x0) != nx
        error("size(x0) must match the number of states of sys")
    elseif size(u(1,x0)) != (nu,) && size(u(1,x0)) != (nu,1)
        error("return value of u must be of size nu")
    end

    dt = Float64(t[2] - t[1])
    if !iscontinuous(sys) || method == :zoh
        if iscontinuous(sys)
            dsys = c2d(sys, dt, :zoh)[1]
        else
            if sys.Ts != dt
                error("Time vector must match sample time for discrete system")
            end
            dsys = sys
        end
    else
        dsys, x0map = c2d(sys, dt, :foh)
        x0 = x0map*[x0; u(1,x0)]
    end
    x,uout = ltitr(dsys.A, dsys.B, u, length(t), map(Float64,x0))
    y = (sys.C*(x.') + sys.D*(uout.')).'
    return y, t, x, uout
end

lsim(sys::TransferFunction, u, t, args...) = lsim(ss(sys), u, t, args...)

function lsim(sys::TransferFunction{SisoGeneralized}, u, t)
    ny, nu = size(sys)
    if !any(size(u) .== [(length(t), nu) (length(t),)])
        error("u must be of size (length(t), nu)")
    end
    ny, nu = size(sys)
    y = Array{Float64}(length(t),ny)
    for i = 1:nu
        for o = 1:ny
            dt = Float64(t[2]-t[1])
            y[:,o] += lsimabstract(sys.matrix[o,i], u[:,i], dt, t[end])
        end
    end
    return y, t
end

@doc """`ltitr(A, B, u[,x0])`

`ltitr(A, B, u::Function, iters[,x0])`

Simulate the discrete time system `x[k + 1] = A x[k] + B u[k]`, returning `x`.
If `x0` is not provided, a zero-vector is used.

If `u` is a function, then `u(i,x)` is called to calculate the control signal every iteration. This can be used to provide a control law such as state feedback `u=-Lx` calculated by `lqr`. In this case, an integrer `iters` must be provided that indicates the number of iterations.
""" ->
function ltitr{T}(A::Matrix{T}, B::Matrix{T}, u::AbstractVecOrMat{T},
        x0::VecOrMat{T}=zeros(T, size(A, 1), 1))
    n = size(u, 1)
    x = Array(T, size(A, 1), n)
    for i=1:n
        x[:,i] = x0
        x0 = A * x0 + B * u[i,:].'
    end
    return x.'
end


function ltitr{T}(A::Matrix{T}, B::Matrix{T}, u::Function, iters::Int,
    x0::VecOrMat{T}=zeros(T, size(A, 1), 1))
    x = Array(T, size(A, 1), iters)
    uout = Array(T, size(B, 2), iters)

    for i=1:iters
        x[:,i] = x0
        uout[:,i] = u(i,x0)
        x0 = A * x0 + B * uout[:,i]
    end
    return x.', uout.'
end

# HELPERS:

# TODO: This is a poor heuristic to estimate a "good" time vector to use for
# simulation, in cases when one isn't provided.
function _default_time_vector(sys::LTISystem, Tf::Real=-1)
    Ts = _default_Ts(sys)
    if Tf == -1
        Tf = 100*Ts
    end
    return linspace(0, Tf, round(Int, Tf/Ts))
end

function _default_Ts(sys::LTISystem)
    if !iscontinuous(sys)
        Ts = sys.Ts
    elseif !isstable(sys)
        Ts = 0.05
    else
        ps = pole(sys)
        r = minimum([abs(real(ps));0])
        if r == 0.0
            r = 1.0
        end
        Ts = 0.07/r
    end
    return Ts
end

_default_Ts(sys::TransferFunction{SisoGeneralized}) = 0.07

#TODO a reasonable check
_issmooth(u::Function) = false

# Determine if a signal is "smooth"
function _issmooth(u, thresh::AbstractFloat=0.75)
    u = [zeros(1, size(u, 2)); u]       # Start from 0 signal always
    dist = maximum(u) - minimum(u)
    du = abs(diff(u))
    return !isempty(du) && all(maximum(du) <= thresh*dist)
end
