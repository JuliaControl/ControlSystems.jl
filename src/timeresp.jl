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

Genereate a workspace object for use with the in-place function [`lsim!`](@ref).
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

Calculate the step response of system `sys`. If the final time `tfinal` or time
vector `t` is not provided, one is calculated based on the system pole
locations. The return value is a structure of type `SimResult` that can be plotted or destructured as `y, t, x = result`.

`y` has size `(ny, length(t), nu)`, `x` has size `(nx, length(t), nu)`"""
function Base.step(sys::AbstractStateSpace, t::AbstractVector; method=:cont, kwargs...)
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

Calculate the impulse response of system `sys`. If the final time `tfinal` or time
vector `t` is not provided, one is calculated based on the system pole
locations. The return value is a structure of type `SimResult` that can be plotted or destructured as `y, t, x = result`.

`y` has size `(ny, length(t), nu)`, `x` has size `(nx, length(t), nu)`"""
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

Calculate the time response of system `sys` to input `u`. If `x0` is ommitted,
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

Continuous-time systems are simulated using an ODE solver if `u` is a function. If `u` is an array, the system is discretized (with `method=:zoh` by default) before simulation. For a lower-level inteface, see `?Simulator` and `?solve`

`u` can be a function or a matrix/vector of precalculated control signals.
If `u` is a function, then `u(x,i)` (`u(x,t)`) is called to calculate the control signal every iteration (time instance used by solver). This can be used to provide a control law such as state feedback `u(x,t) = -L*x` calculated by `lqr`.
To simulate a unit step at `t=t₀`, use `(x,i)-> Ts*i ≥ t₀`, for a ramp, use `(x,i)-> i*Ts`, for a step at `t=5`, use (x,i)-> (i*Ts >= 5) etc.

For maximum performance, see function [`lsim!`](@ref), avaialable for discrete-time systems only.

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

u(x,t) = -L*x # Form control law
t  = 0:0.1:5
x0 = [1,0]
y, t, x, uout = lsim(sys,u,t,x0=x0)
plot(t,x', lab=["Position" "Velocity"], xlabel="Time [s]")
```
"""
function lsim(sys::AbstractStateSpace, u::AbstractVecOrMat, t::AbstractVector;
        x0::AbstractVecOrMat=zeros(Bool, nstates(sys)), method::Symbol=:zoh)
    ny, nu = size(sys)
    nx = sys.nx

    if length(x0) != nx
        error("size(x0) must match the number of states of sys")
    end
    if size(u) != (nu, length(t))
        error("u must be of size (nu, length(t))")
    end

    dt = Float64(t[2] - t[1])
    if !all(x -> x ≈ dt, diff(t))
        error("time vector t must be uniformly spaced")
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
        if sys.Ts != dt
            error("Time vector must match sample time of discrete-time system")
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

function lsim(sys::AbstractStateSpace, u::Function, t::AbstractVector;
        x0::AbstractVecOrMat=zeros(Bool, nstates(sys)), method::Symbol=:cont, alg = Tsit5(), kwargs...)
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
            simsys = c2d(sys, dt, :zoh)
        else
            if sys.Ts != dt
                error("Time vector must match sample time for discrete system")
            end
            simsys = sys
        end
        x,uout = ltitr(simsys.A, simsys.B, u, t, T.(x0))
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


lsim(sys::TransferFunction, args...; kwargs...) = lsim(ss(sys), args...; kwargs...)


"""
    ltitr(A, B, u[,x0])
    ltitr(A, B, u::Function, iters[,x0])

Simulate the discrete time system `x[k + 1] = A x[k] + B u[k]`, returning `x`.
If `x0` is not provided, a zero-vector is used.

The type of `x0` determines the matrix structure of the returned result,
e.g, `x0` should prefereably not be a sparse vector.

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
    size(x, 2) == length(t) || throw(ArgumentError("Inconsitent lengths of workspace cache and t"))
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
    dt = _default_dt(sys)
    if tfinal == -1
        tfinal = 200dt
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
