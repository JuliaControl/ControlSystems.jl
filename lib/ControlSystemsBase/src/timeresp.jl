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
one is calculated based on the system pole locations. 

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
Base.step(sys::TransferFunction, t::AbstractVector; kwargs...) = step(ss(sys), t::AbstractVector; kwargs...)

"""
    y, t, x = impulse(sys[, tfinal])
    y, t, x = impulse(sys[, t])

Calculate the response of the system `sys` to an impulse at time `t = 0`. 
For continous-time systems, the impulse is a unit Dirac impulse. 
For discrete-time systems, the impulse lasts one sample and has magnitude `1/Ts`. 
If the final time `tfinal` or time vector `t` is not provided, 
one is calculated based on the system pole locations. 

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
impulse(sys::TransferFunction, t::AbstractVector; kwargs...) = impulse(ss(sys), t; kwargs...)

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

    if !iscontinuous(sys) || method === :zoh
        if iscontinuous(sys)
            simsys = c2d(sys, dt, :zoh)
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


lsim(sys::TransferFunction, args...; kwargs...) = lsim(ss(sys), args...; kwargs...)


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

# TODO: This is a poor heuristic to estimate a "good" time vector to use for
# simulation, in cases when one isn't provided.
function _default_time_vector(sys::LTISystem, tfinal::Real=-1)
    dt = _default_dt(sys) # This is set small enough to resolve the fastest dynamics. 
    if tfinal == -1
        # The final simulation time should be chosen to also let the slowest dynamics run its course. 
        # Integrating systems pose a problem for this heuristic, and we sort out any poles that appear to be integrators.
        if hasmethod(poles, typeof((sys,)))
            # DelaySystem does not define poles
            ws = abs.(poles(sys))
            ω0_min = minimum(w for w in ws if w > 1e-6; init=dt)
            dt_slow = round(1/(2ω0_min), sigdigits=2)
            tfinal = max(200dt, dt_slow)
            tfinal = min(tfinal, 100_000*dt)
        else
            tfinal = 200dt
        end
    elseif iscontinuous(sys)
        dt = min(dt, tfinal/200)
    end
    return 0:dt:tfinal
end

function _default_dt(sys::LTISystem)
    if isdiscrete(sys)
        return sys.Ts
    elseif all(iszero, poles(sys)) # Static or pure integrators
        return 0.05
    else
        ω0_max = maximum(abs.(poles(sys)))
        dt = round(1/(12*ω0_max), sigdigits=2)
        return dt
    end
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
    settlingtimeind == length(res.t) && @warn "System might not have settled within the simulation time"
    settlingtime = res.t[settlingtimeind] + Ts
    op = direction == 1 ? (>) : (<)
    i10 = findfirst(op.(y, y0 + risetime_th[1] * stepsize * direction))
    i90 = findfirst(op.(y, y0 + risetime_th[2] * stepsize * direction))
    risetime = res.t[i90] - res.t[i10]
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

