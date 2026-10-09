# Functions for calculating time response of a system

# XXX : `step` is a function in Base, with a different meaning than it has
# here. This shouldn't be an issue, but it might be.
struct LsimWorkspace{T}
    x::Matrix{T}
    u::Matrix{T}
    y::Matrix{T}
end
function LsimWorkspace{T}(ny::Int, nu::Int, nx::Int, N::Int) where T
    x = Matrix{T}(undef, nx, N)
    u = Matrix{T}(undef, nu, N)
    y = Matrix{T}(undef, ny, N)
    LsimWorkspace{T}(x, u, y)
end

"""
    LsimWorkspace(sys::AbstractStateSpace, N::Int)
    LsimWorkspace(sys::AbstractStateSpace, u::AbstractMatrix)
    LsimWorkspace{T}(ny, nu, nx, N)

Generate a workspace object for use with the in-place function [`lsim!`](@ref).
`sys` is the discrete-time system to be simulated and `N` is the number of time steps, alternatively, the input `u` can be provided instead of `N`.
Note: for threaded applications, create one workspace object per thread. 
"""
function LsimWorkspace(sys::AbstractStateSpace, N::Int)
    T = numeric_type(sys)
    x = Matrix{T}(undef, sys.nx, N)
    u = Matrix{T}(undef, sys.nu, N)
    y = Matrix{T}(undef, sys.ny, N)
    LsimWorkspace(x, u, y)
end

LsimWorkspace(sys::AbstractStateSpace, u::AbstractMatrix) = LsimWorkspace(sys, size(u, 2))

"""
    y, t, x = step(sys[, tfinal])
    y, t, x = step(sys[, t])

Calculate the response of the system `sys` to a unit step at time `t = 0`. 
If the final time `tfinal` or time vector `t` is not provided, 
one is calculated based on the system pole locations: the sample interval
resolves the fastest pole and the final time is chosen such that the slowest
stable mode has decayed, or, for an unstable system, such that the fastest
growing mode has increased by a factor of about ``e^5``. Integrators are
excluded from this computation.

The return value is a structure of type `SimResult`. 
A `SimResul` can be plotted by `plot(result)`, 
or destructured as `y, t, x = result`. 

`y` has size `(ny, length(t), nu)`, `x` has size `(nx, length(t), nu)`

See also [`stepinfo`](@ref) and [`lsim`](@ref).
"""
function Base.step(sys::AbstractStateSpace, t::AbstractVector; method=:zoh, kwargs...)
    T = promote_type(eltype(sys.A), Float64)
    ny, nu = size(sys)
    nx = nstates(sys)
    u = let u_element = [one(eltype(t))] # to avoid allocating this multiple times
        (x,t)->u_element
    end
    x0 = zeros(T, nx)
    if nu == 1
        y, tout, x, uout = lsim(sys, u, t; x0, method, kwargs...)
    else
        x = Array{T}(undef, nx, length(t), nu)
        y = Array{T}(undef, ny, length(t), nu)
        for i=1:nu
            y[:,:,i], tout, x[:,:,i], uout = lsim(sys[:,i], u, t; x0, method, kwargs...)
        end
    end
    return SimResult(y, t, x, uout, sys)
end

Base.step(sys::LTISystem, tfinal::Real; kwargs...) = step(sys, _default_time_vector(sys, tfinal); kwargs...)
Base.step(sys::LTISystem; kwargs...) = step(sys, _default_time_vector(sys); kwargs...)
Base.step(sys::TransferFunction, t::AbstractVector; kwargs...) = step(ss(sys, minimal=numeric_type(sys) isa BlasFloat), t::AbstractVector; kwargs...)

"""
    y, t, x = impulse(sys[, tfinal])
    y, t, x = impulse(sys[, t])

Calculate the response of the system `sys` to an impulse at time `t = 0`. 
For continous-time systems, the impulse is a unit Dirac impulse. 
For discrete-time systems, the impulse lasts one sample and has magnitude `1/Ts`. 
If the final time `tfinal` or time vector `t` is not provided, 
one is calculated based on the system pole locations: the sample interval
resolves the fastest pole and the final time is chosen such that the slowest
stable mode has decayed, or, for an unstable system, such that the fastest
growing mode has increased by a factor of about ``e^5``. Integrators are
excluded from this computation.

The return value is a structure of type `SimResult`. 
A `SimResul` can be plotted by `plot(result)`, 
or destructured as `y, t, x = result`.

`y` has size `(ny, length(t), nu)`, `x` has size `(nx, length(t), nu)`

See also [`lsim`](@ref).
"""
function impulse(sys::AbstractStateSpace, t::AbstractVector; kwargs...)
    T = promote_type(eltype(sys.A), Float64)
    ny, nu = size(sys)
    nx = nstates(sys)
    if iscontinuous(sys) #&& method === :cont
        u = (x,t) -> [zero(T)]
        # impulse response equivalent to unforced response with x0 = B.
        imp_sys = ss(sys.A, zeros(T, nx, nu), sys.C, 0)
        x0s = sys.B
    else
        u_element = [zero(T)]
        u = (x,i) -> (i == t[1] ? [one(T)]/sys.Ts : u_element)
        x0s = zeros(T, nx, nu)
    end
    if nu == 1 # Why two cases # QUESTION: Not type stable?
        y, t, x, uout = lsim(sys, u, t; x0=x0s[:], kwargs...)
    else
        x = Array{T}(undef, nx, length(t), nu)
        y = Array{T}(undef, ny, length(t), nu)
        for i=1:nu
            y[:,:,i], t, x[:,:,i], uout = lsim(sys[:,i], u, t; x0=x0s[:,i], kwargs...)
        end
    end
    return SimResult(y, t, x, uout, sys)
end

impulse(sys::LTISystem, tfinal::Real; kwargs...) = impulse(sys, _default_time_vector(sys, tfinal); kwargs...)
impulse(sys::LTISystem; kwargs...) = impulse(sys, _default_time_vector(sys); kwargs...)
impulse(sys::TransferFunction, t::AbstractVector; kwargs...) = impulse(ss(sys, minimal=numeric_type(sys) isa BlasFloat), t; kwargs...)

"""
    result = lsim(sys, u[, t]; x0, method])
    result = lsim(sys, u::Function, t; x0, method)

Calculate the time response of system `sys` to input `u`. If `x0` is omitted,
a zero vector is used.

The result structure contains the fields `y, t, x, u` and can be destructured automatically by iteration, e.g.,
```julia
y, t, x, u = result
```
`result::SimResult` can also be plotted directly:
```julia
plot(result, plotu=true, plotx=false)
```
`y`, `x`, `u` have time in the second dimension. Initial state `x0` defaults to zero.

Continuous-time systems are simulated using an ODE solver if `u` is a function (requires using ControlSystems). If `u` is an array, the system is discretized (with `method=:zoh` by default) before simulation. For a lower-level interface, see `?Simulator` and `?solve`. For continuous-time systems, keyword arguments are forwarded to the ODE solver. By default, the option `dtmax = t[2]-t[1]` is used to prevent the solver from stepping over discontinuities in `u(x, t)`. This prevents the solver from taking too large steps, but may also slow down the simulation when `u` is smooth. To disable this behavior, set `dtmax = Inf`.

`u` can be a function or a *matrix* of precalculated control signals and must have dimensions `(nu, length(t))`.
If `u` is a function, then `u(x,i)` (for discrete systems) or `u(x,t)` (for continuous ones) is called to calculate the control signal at every iteration (time instance used by solver). This can be used to provide a control law such as state feedback `u(x,t) = -L*x` calculated by `lqr`.
To simulate a unit step at `t=t₀`, use `(x,t)-> t ≥ t₀`, for a ramp, use `(x,t)-> t`, for a step at `t=5`, use `(x,t)-> (t >= 5)` etc.

*Note:* The function `u` will be called once before simulating to verify that it returns an array of the correct dimensions. This can cause problems if `u` is stateful or has other side effects. You can disable this check by passing `check_u = false`.

For maximum performance, see function [`lsim!`](@ref), available for discrete-time systems only.

Usage example:
```julia
using ControlSystems
using LinearAlgebra: I
using Plots

A = [0 1; 0 0]
B = [0;1]
C = [1 0]
sys = ss(A,B,C,0)
Q = I
R = I
L = lqr(sys,Q,R)

u(x,t) = -L*x # Form control law
t  = 0:0.1:5
x0 = [1,0]
y, t, x, uout = lsim(sys,u,t,x0=x0)
plot(t,x', lab=["Position" "Velocity"], xlabel="Time [s]")

# Alternative way of plotting
res = lsim(sys,u,t,x0=x0)
plot(res)
```
"""
function lsim(sys::AbstractStateSpace, u::AbstractVecOrMat, t::AbstractVector;
        x0::AbstractVecOrMat=zeros(eltype(u), nstates(sys)), method::Symbol=:zoh)
    ny, nu = size(sys)
    nx = sys.nx
    T = typeof(u)

    length(x0) == nx ||
        error("x0 must have length $nx: got length $(length(x0))")

    if size(u) != (nu, length(t))
        if (size(u) == (length(t), nu)) || (u isa AbstractVector && length(u) == length(t))
            @warn("u should be a matrix of size ($nu, $(length(t))): got an u of type $(typeof(u)) and size $(size(u)). An attempt at using the transpose of u will be performed (this may fail if it can not be done in a type stable way, in which case you get a TypeError). To silence this warning, use the correct input dimension.")
            u = copy(transpose(u))::T
        else
            error("u must have size ($nu, $(length(t))): got size $(size(u))")
        end
    end

    dt = Float64(t[2] - t[1])
    if !all(x -> x ≈ dt, diff(t))
        error("Time vector t must be uniformly spaced")
    end

    # Handle pure D system (no states)
    if nx == 0
        x = Matrix{eltype(u)}(undef, 0, length(t)) # no states
        y = sys.D * u
        dsys = sys
        return SimResult(y, t, x, u, dsys)
    end

    if iscontinuous(sys)
        if method === :zoh
            dsys = c2d(sys, dt, :zoh)
        elseif method === :foh
            dsys, x0map = c2d_x0map(sys, dt, :foh)
            x0 = x0map*[x0; u[:,1]]
        else
            error("Unsupported discretization method: $method")
        end
    else
        if !(sys.Ts ≈ dt)
            error("Time vector must match the sample time of the discrete-time system, $(sys.Ts): got $dt")
        end
        dsys = sys
    end

    x = ltitr(dsys.A, dsys.B, u, x0)
    y = sys.C*x
    if !iszero(sys.D)
        mul!(y, sys.D, u, 1, 1)
    end
    return SimResult(y, t, x, u, dsys) # saves the system that actually produced the simulation
end

function lsim(sys::AbstractStateSpace{<:Discrete}, u::AbstractVecOrMat; kwargs...)
    nu = sys.nu
    if size(u, 1) != nu
        if u isa AbstractVector && sys.nu == 1
            # The isa Array is a safeguard against type instability due to the copy(u')
            @warn("u should be a row-vector of size (1, $(length(u))): got a regular vector of size $(size(u)). The transpose of u will be used. To silence this warning, use the correct input dimension.")
            u = copy(transpose(u))
        else
            error("u must be a matrix of size (nu, length(t)) where the number of inputs nu=$nu, got size u = $(size(u))")
        end
    end
    t = range(0, length=size(u, 2), step=sys.Ts)
    lsim(sys, u, t; kwargs...)
end

@deprecate lsim(sys, u, t, x0) lsim(sys, u, t; x0)
@deprecate lsim(sys, u, t, x0, method) lsim(sys, u, t; x0, method)

function lsim(sys::AbstractStateSpace, u::Function, tfinal::Real; kwargs...)
    t = _default_time_vector(sys, tfinal)
    lsim(sys, u, t; kwargs...)
end

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

# This method is less specific than the lsim in ControlSystems that specifies Continuous timeevol for sys, hence, if ControlSystems is loaded, ControlSystems.lsim will take precedence over this
function lsim(sys::AbstractStateSpace, u::Function, t::AbstractVector;
        x0::AbstractVecOrMat=zeros(Bool, nstates(sys)), method::Symbol=:cont, check_u = true, kwargs...)
    ny, nu = size(sys)
    nx = sys.nx
    if length(x0) != nx
        error("x0 must have length $nx: got length $(length(x0))")
    end
    if check_u
        u0 = u(x0,t[1])
        if !(u0 isa Number && nu == 1) && (size(u0) != (nu,) && size(u0) != (nu,1))
            error("u returned by the input function must have size ($nu,): got size(u0) = $(size(u0))")
        end
    end
    T = promote_type(Float64, eltype(x0), numeric_type(sys))

    dt = t[2] - t[1]

    # Handle pure D system (no states)
    if nx == 0
        uout = Matrix{T}(undef, nu, length(t))
        for i = 1:length(t)
            uout[:, i] = u(T[], t[i])
        end
        x = Matrix{T}(undef, 0, length(t)) # no states
        y = sys.D * uout
        simsys = sys
        return SimResult(y, t, x, uout, simsys)
    end

    if !iscontinuous(sys) || method ∈ (:zoh, :tustin, :foh, :fwdeuler)
        if iscontinuous(sys)
            simsys = c2d(sys, dt, method)
        else
            if !(sys.Ts ≈ dt)
                error("Time vector interval ($dt) must match sample time for discrete system ($(sys.Ts))")
            end
            simsys = sys
        end
        x,uout = ltitr(simsys.A, simsys.B, u, t, T.(x0))
    else
        throw(MethodError(lsim, (sys, u, t)))
    end
    y = sys.C*x
    if !iszero(sys.D)
        mul!(y, sys.D, uout, 1, 1)
    end
    return SimResult(y, t, x, uout, simsys) # saves the system that actually produced the simulation
end


lsim(sys::TransferFunction, args...; kwargs...) = lsim(ss(sys, minimal=numeric_type(sys) isa BlasFloat), args...; kwargs...)


"""
    ltitr(A, B, u[,x0])
    ltitr(A, B, u::Function, iters[,x0])

Simulate the discrete time system `x[k + 1] = A x[k] + B u[k]`, returning `x`.
If `x0` is not provided, a zero-vector is used.

The type of `x0` determines the matrix structure of the returned result,
e.g, `x0` should preferably not be a sparse vector.

If `u` is a function, then `u(x,i)` is called to calculate the control signal every iteration. This can be used to provide a control law such as state feedback `u=-Lx` calculated by `lqr`. In this case, an integrer `iters` must be provided that indicates the number of iterations.
"""
function ltitr(A::AbstractMatrix, B::AbstractMatrix, u::AbstractVecOrMat,
    x0::AbstractVecOrMat=zeros(eltype(A), size(A, 1)))

    T = promote_type(LinearAlgebra.promote_op(LinearAlgebra.matprod, eltype(A), eltype(x0)),
                  LinearAlgebra.promote_op(LinearAlgebra.matprod, eltype(B), eltype(u)))

    n = size(u, 2)
    # Using similar instead of Matrix{T} to allow for CuArrays to be used.
    # This approach is problematic if x0 is sparse for example, but was considered
    # to be good enough for now
    x = similar(x0, T, (length(x0), n))
    ltitr!(x, A, B, u, x0)
end

function ltitr(A::AbstractMatrix{T}, B::AbstractMatrix{T}, u::Function, t,
    x0::AbstractVecOrMat=zeros(T, size(A, 1))) where T
    iters = length(t)
    x = similar(A, size(A, 1), iters)
    uout = similar(A, size(B, 2), iters)
    ltitr!(x, uout, A, B, u, t, x0)
end

## In-place
@views function ltitr!(x, A::AbstractMatrix, B::AbstractMatrix, u::AbstractVecOrMat,
    x0::AbstractVecOrMat=zeros(eltype(A), size(A, 1)))
    x[:,1] .= x0
    n = size(u, 2)
    mul!(x[:, 2:end], B, u[:, 1:end-1]) # Do all multiplications B*u[:,k] to save view allocations

    for k=1:n-1
        mul!(x[:, k+1], A, x[:,k], true, true)
    end
    return x
end

@views function ltitr!(x, uout, A::AbstractMatrix{T}, B::AbstractMatrix{T}, u::Function, t,
    x0::AbstractVecOrMat=zeros(T, size(A, 1))) where T
    size(x,2) == length(t) == size(uout, 2) || throw(ArgumentError("Inconsistent array sizes."))
    x[:, 1] .= x0
    @inbounds for i=1:length(t)-1
        xi = x[:,i]
        xp = x[:,i+1]
        uout[:,i] .= u(xi,t[i])
        mul!(xp, B, uout[:,i])
        mul!(xp, A, xi, true, true)
    end
    uout[:,end] .= u(x[:,end],t[end])
    return x, uout
end

function lsim!(ws::LsimWorkspace, sys::AbstractStateSpace{<:Discrete}, u::AbstractVecOrMat; kwargs...)
    t = range(0, length=size(u, 2), step=sys.Ts)
    lsim!(ws, sys, u, t; kwargs...)
end

"""
    res = lsim!(ws::LsimWorkspace, sys::AbstractStateSpace{<:Discrete}, u, [t]; x0)

In-place version of [`lsim`](@ref) that takes a workspace object created by calling [`LsimWorkspace`](@ref).
*Notice*, if `u` is a function, `res.u === ws.u`. If `u` is an array, `res.u === u`.
"""
function lsim!(ws::LsimWorkspace{T}, sys::AbstractStateSpace{<:Discrete}, u, t::AbstractVector;
        x0::AbstractVecOrMat=zeros(eltype(T), nstates(sys))) where T

    x, y = ws.x, ws.y
    size(x, 2) == length(t) || throw(ArgumentError("Inconsistent lengths of workspace cache and t"))
    arr = u isa AbstractArray
    if arr
        copyto!(ws.u, u)
        ltitr!(x, sys.A, sys.B, u, x0)
    else # u is a function
        ltitr!(x, ws.u, sys.A, sys.B, u, t, x0)
    end
    mul!(y, sys.C, x)
    if !iszero(sys.D)
        mul!(y, sys.D, arr ? u : ws.u, 1, 1)
    end
    if arr
        SimResult(y, t, x, u, sys)
    else
        SimResult(y, t, x, ws.u, sys)
    end
end

# HELPERS:

"""
    t = _default_time_vector(sys, tfinal = -1)

Compute a time vector for simulation of `sys` when none is provided.

The sample interval `dt` resolves the fastest pole, see `_default_dt`. If `tfinal` is not provided, it is computed from the poles by `_default_tfinal`, such that the slowest stable mode has decayed. The number of samples is limited to 100 001: for a continuous-time system, `dt` is increased if required, for a discrete-time system, whose sample interval is fixed, `tfinal` is reduced and a warning is emitted.

The time vector is invariant under `time_scale` up to the rounding of `dt` and `tfinal` to two significant digits.
"""
function _default_time_vector(sys::LTISystem, tfinal::Real=-1)
    dt = _default_dt(sys) # This is set small enough to resolve the fastest dynamics.
    if tfinal == -1
        if hasmethod(poles, typeof((sys,)))
            # DelaySystem does not define poles
            tfinal = _default_tfinal(sys, dt)
            nmax = 100_000
            if tfinal > nmax*dt
                if isdiscrete(sys)
                    @warn "The default final time $tfinal of the simulation of the discrete-time system requires more than $nmax samples, it is reduced to $(nmax*dt). Provide the final time or a time vector to simulate for a longer time."
                    tfinal = nmax*dt
                else
                    dt = round(tfinal/nmax, RoundUp, sigdigits=2)
                end
            end
        else
            tfinal = 200dt
        end
    elseif iscontinuous(sys)
        dt = min(dt, tfinal/200)
    end
    return 0:dt:tfinal
end

"""
    p, tol = _nonintegrator_poles(sys)

Return the poles of `sys` that are not integrators, mapped to the s-plane by s = log(z)/Ts for discrete-time systems, together with a tolerance `tol` below which the real part of such a pole is considered zero. A pole at z = 0 is mapped to s = -∞.

Poles in the origin (z = 1 in discrete time) are classified as integrators with the tolerance of `count_eigval_multiplicity`, which accounts for the larger rounding errors of multiple poles. The rounding errors of the eigenvalues of `A` are proportional to the norm of `A`, which may be much larger than the magnitude of the computed poles, e.g., if all poles are in the origin. The norm of `A` is therefore included in the scale of the tolerance.
"""
function _nonintegrator_poles(sys::LTISystem)
    p = poles(sys)
    location = iscontinuous(sys) ? 0 : 1
    scale = float(maximum(abs, p, init=0.0))
    if sys isa AbstractStateSpace && !isempty(sys.A)
        scale = max(scale, opnorm(sys.A, 1))
    end
    nint, tol = count_eigval_multiplicity(p, location; scale)
    if nint > 0
        p = filter(p -> abs(p - location) > tol, p)
    end
    tol_axis = 100*eps(float(real(eltype(p))))*scale # Tolerance for a simple pole on the imaginary axis (unit circle)
    if isdiscrete(sys)
        return log.(complex.(p)) ./ sys.Ts, tol_axis/sys.Ts
    end
    p, tol_axis
end

"""
    tfinal = _default_tfinal(sys, dt)

Final time of the default simulation time vector with sample interval `dt`, computed from the poles that are not integrators:
- If there are unstable poles, `tfinal = 5τ` where `τ` is the time constant of the fastest growing mode, such that the growth is visible but the response does not overflow.
- Otherwise, `tfinal = μ + 6s` with `μ = Σ τᵢ` and `s = √(Σ τᵢ²)`, where `τᵢ` are the time constants of the stable poles. The impulse response of a series connection of first-order systems with time constants `τᵢ` equals the probability density of a sum of independent exponentially distributed random variables with means `τᵢ`, whose mean and standard deviation are `μ` and `s`. For a single pole, `tfinal = 7τ`, after which a fraction `exp(-7) < 1e-3` of the initial deviation remains. For poles of high multiplicity, `s` accounts for the polynomial factors of the response. The time constant of a complex pole is the inverse of the magnitude of its real part, i.e., the time constant of the envelope of the oscillation. Each pole of a complex-conjugate pair is counted, which yields `tfinal ≈ 10.5τ` for a single pair. Counting each pair once would yield a final time that depends discontinuously on the poles, since rounding errors split a multiple real pole into complex-conjugate pairs. The time constants are bounded from below by `dt`, since a discrete-time mode does not decay faster than in one sample interval.
- `tfinal` is at least five periods of the slowest undamped oscillatory mode.
- If all poles are integrators, or the system has no poles, the system has no characteristic time and `tfinal = 200dt`.
"""
function _default_tfinal(sys::LTISystem, dt)
    p, tol = _nonintegrator_poles(sys)
    isempty(p) && return 200dt
    σ = -real.(p) # Decay rates
    growth_rate = -minimum(σ)
    if growth_rate > tol
        return round(5*max(1/growth_rate, dt), sigdigits=2)
    end
    stable = σ .> tol
    τ = [max(1/σ[i], dt) for i in eachindex(p) if stable[i]]
    tfinal = sum(τ, init=0.0) + 6*sqrt(sum(abs2, τ, init=0.0))
    ω_undamped = minimum((abs(p[i]) for i in eachindex(p) if !stable[i]), init=Inf)
    if isfinite(ω_undamped)
        tfinal = max(tfinal, 5*2π/ω_undamped)
    end
    round(tfinal, sigdigits=2)
end

"""
    dt = _default_dt(sys)

Default sample interval for simulation of `sys`. For a discrete-time system, `dt = sys.Ts`. For a continuous-time system, `dt = 1/(12 max|pᵢ|)` resolves the fastest pole, where integrators are excluded. If all poles are integrators, or the system has no poles, the system has no characteristic time and `dt = 0.05`.
"""
function _default_dt(sys::LTISystem)
    isdiscrete(sys) && return sys.Ts
    p, _ = _nonintegrator_poles(sys)
    isempty(p) && return 0.05
    round(1/(12*maximum(abs, p)), sigdigits=2)
end

"""
    StepInfo

Computed using [`stepinfo`](@ref)

# Fields:
- `y0`: The initial value of the step response.
- `yf`: The final value of the step response.
- `stepsize`: The size of the step.
- `peak`: The peak value of the step response.
- `peaktime`: The time at which the peak occurs.
- `overshoot`: The overshoot of the step response.
- `settlingtime`: The time at which the step response has settled to within `settling_th` of the final value.
- `settlingtimeind::Int`: The index at which the step response has settled to within `settling_th` of the final value.
- `risetime`: The time at which the response rises from `risetime_th[1]` to `risetime_th[2]` of the final value 
- `i10::Int`: The index at which the response reaches `risetime_th[1]`
- `i90::Int`: The index at which the response reaches `risetime_th[2]`
- `res::SimResult{SR}`: The simulation result used to compute the step response characteristics.
- `settling_th`: The threshold used to compute `settlingtime` and `settlingtimeind`.
- `risetime_th`: The thresholds used to compute `risetime`, `i10`, and `i90`.
"""
struct StepInfo{SR}
    y0
    yf
    stepsize
    peak
    peaktime
    overshoot
    lowerpeak
    lowerpeakind
    undershoot
    settlingtime
    settlingtimeind::Int
    risetime
    i10::Int
    i90::Int
    res::SimResult{SR}
    settling_th
    risetime_th
end

"""
    stepinfo(res::SimResult; y0 = nothing, yf = nothing, settling_th = 0.02, risetime_th = (0.1, 0.9))

Compute the step response characteristics for a simulation result. The following information is computed and stored in a [`StepInfo`](@ref) struct:
- `y0`: The initial value of the response
- `yf`: The final value of the response
- `stepsize`: The size of the step
- `peak`: The peak value of the response
- `peaktime`: The time at which the peak occurs
- `overshoot`: The percentage overshoot of the response
- `undershoot`: The percentage undershoot of the response. If the step response never reaches below the initial value, the undershoot is zero.
- `settlingtime`: The time at which the response settles within `settling_th` of the final value
- `settlingtimeind`: The index at which the response settles within `settling_th` of the final value
- `risetime`: The time at which the response rises from `risetime_th[1]` to `risetime_th[2]` of the final value


# Arguments:
- `res`: The result from a simulation using [`step`](@ref) (or [`lsim`](@ref))
- `y0`: The initial value, if not provided, the first value of the response is used.
- `yf`: The final value, if not provided, the last value of the response is used. The simulation must have reached steady-state for an automatically computed value to make sense. If the simulation has not reached steady state, you may provide the final value manually.
- `settling_th`: The threshold for computing the settling time. The settling time is the time at which the response settles within `settling_th` of the final value.
- `risetime_th`: The lower and upper threshold for computing the rise time. The rise time is the time at which the response rises from `risetime_th[1]` to `risetime_th[2]` of the final value.

# Example:
```julia
G = tf([1], [1, 1, 1])
res = step(G, 15)
si = stepinfo(res)
plot(si)
```
"""
function stepinfo(res::SimResult; y0 = nothing, yf = nothing, settling_th = 0.02, risetime_th = (0.1, 0.9))
    issiso(res) || throw(ArgumentError("stepinfo only supports SISO systems"))
    y = res.y[1, :]
    y0 === nothing && (y0 = y[1])
    yf === nothing && (yf = y[end])
    Ts = res.t[2] - res.t[1]
    direction = sign(yf - y0)
    stepsize = abs(yf - y0)
    peak, peakind = direction == 1 ? findmax(y) : findmin(y)
    lowerpeak, lowerpeakind = direction == 1 ? findmin(y) : findmax(y)
    peaktime = res.t[peakind]
    overshoot = direction * 100 * (peak - yf) / stepsize
    undershoot = direction * 100 * (y0 - lowerpeak) / stepsize
    undershoot > 0 ? undershoot : zero(undershoot)
    settlingtimeind = findlast(abs(y-yf) > settling_th * stepsize for y in y)
    settlingtimeind === nothing && (settlingtimeind = length(res.t))
    settlingtimeind == length(res.t) && @warn "System might not have settled within the simulation time"
    settlingtime = res.t[settlingtimeind] + Ts
    op = direction == 1 ? (>) : (<)
    i10 = findfirst(op.(y, y0 + risetime_th[1] * stepsize * direction))
    i90 = findfirst(op.(y, y0 + risetime_th[2] * stepsize * direction))
    if i10 === nothing || i90 === nothing
        @warn "Response did not reach the requested risetime threshold(s) within the simulation window"
        i10 = i10 === nothing ? 0 : i10
        i90 = i90 === nothing ? 0 : i90
        risetime = oftype(float(res.t[1]), NaN)
    else
        risetime = res.t[i90] - res.t[i10]
    end
    StepInfo(y0, yf, stepsize, peak, peaktime, overshoot, lowerpeak, lowerpeakind, undershoot, settlingtime, settlingtimeind, risetime, i10, i90, res, settling_th, risetime_th)
end

function Base.show(io::IO, info::StepInfo)
    println(io, "StepInfo:")
    @printf(io, "%-15s %8.3f\n", "Initial value:", info.y0)
    @printf(io, "%-15s %8.3f\n", "Final value:", info.yf)
    @printf(io, "%-15s %8.3f\n", "Step size:", info.stepsize)
    @printf(io, "%-15s %8.3f\n", "Peak:", info.peak)
    @printf(io, "%-15s %8.3f s\n", "Peak time:", info.peaktime)
    @printf(io, "%-15s %8.2f %%\n", "Overshoot:", info.overshoot)
    @printf(io, "%-15s %8.2f %%\n", "Undershoot:", info.undershoot)
    @printf(io, "%-15s %8.3f s\n", "Settling time:", info.settlingtime)
    @printf(io, "%-15s %8.3f s\n", "Rise time:", info.risetime)
end

