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

"""
    A, B = linearize(f, x, u, args...)

Linearize dynamics ``ẋ = f(x, u, args...)`` around operating point ``(x,u,args...)`` using ForwardDiff. `args` can be empty, or contain, e.g., parameters and time `(p, t)` like in the SciML interface.
This function can also be used to linearize an output equation `C, D = linearize(h, x, u, args...)`.
"""
function linearize(f, xi::AbstractVector, ui::AbstractVector, args...)
    A = ForwardDiff.jacobian(x -> f(x, ui, args...), xi)
    B = ForwardDiff.jacobian(u -> f(xi, u, args...), ui)
    A, B
end

function _bounds_and_features(sys::HammersteinWienerSystem, plot::Val)
    _bounds_and_features(sys.P.P, plot)
end

# Again we have to do something for default vectors, more or less a copy from timeresp.jl
function _default_dt(sys::HammersteinWienerSystem)
    _default_dt(sys.P.P)
end

# The implementation of lsim for HammersteinWienerSystem is located in ControlSystems.timeresp

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


