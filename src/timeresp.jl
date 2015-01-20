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
Base.step(sys::StateSpace, Tf::Real) = step(sys, _default_time_vector(sys, Tf))
Base.step(sys::StateSpace) = step(sys, _default_time_vector(sys))
Base.step(sys::LTISystem, args...) = step(ss(sys), args...)


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
impulse(sys::StateSpace, Tf::Real) = impulse(sys, _default_time_vector(sys, Tf))
impulse(sys::StateSpace) = impulse(sys, _default_time_vector(sys))
impulse(sys::LTISystem, args...) = impulse(ss(sys), args...)


@doc """`[y, t, x] = lsim(sys, u, t[, x0, method])`

Calculate the time response of system `sys` to input `u`. If `x0` is ommitted,
a zero vector is used.

Continuous time systems are discretized before simulation. By default, the
method is chosen based on the smoothness of the input signal. Optionally, the
`method` parameter can be specified as either `:zoh` or `:foh`.""" ->
function lsim(sys::StateSpace, u::AbstractVecOrMat, t::AbstractVector,
        x0::VecOrMat=zeros(sys.nx, 1), method::Symbol=_issmooth(u) ? :foh : :zoh)
    ny, nu = size(sys) 
    nx = sys.nx

    if length(x0) != nx 
        error("size(x0) must match the number of states of sys")
    elseif !any(size(u) .== [(length(t), nu) (length(t),)])
        error("u must be of size (length(t), nu)")
    end

    dt = float64(t[2] - t[1])
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
    x = ltitr(dsys.A, dsys.B, float64(u), float64(x0))
    y = (sys.C*(x.') + sys.D*(u.')).'
    return y, t, x
end
lsim(sys::LTISystem, args...) = lsim(ss(sys), args...)


@doc """`ltitr(A, B, u[, x0])`

Simulate the discrete time system `x[k + 1] = A x[k] + B u[k]`, returning `x`.
If `x0` is not provided, a zero-vector is used.""" ->
function ltitr{T}(A::Matrix{T}, B::Matrix{T}, u::AbstractVecOrMat{T},
        x0::VecOrMat{T})
    n = size(u, 1)
    x = Array(T, size(A, 1), n)
    for i=1:n
        x[:,i] = x0
        x0 = A * x0 + B * u[i,:].'
    end
    return x.'
end 
ltitr{T}(A::Matrix{T}, B::Matrix{T}, u::AbstractVecOrMat{T}) = 
        ltitr(A, B, u, zeros(T, size(A, 1), 1))


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
        r = minimum(abs(real(ps)))
        if r == 0.0
            r = 1.0
        end
        Ts = 0.07/r
    end
    return Ts
end

# Determine if a signal is "smooth"
function _issmooth(u::Array, thresh::FloatingPoint=0.75)
    u = [zeros(1, size(u, 2)); u]       # Start from 0 signal always
    range = maximum(u) - minimum(u)
    du = abs(diff(u))
    return !isempty(du) && all(maximum(du) <= thresh*range)
end
