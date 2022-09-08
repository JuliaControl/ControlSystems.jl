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

```
           y▲   ────── upper
            │  /
            │ /
            │/
  ──────────┼────────► u
           /│   
          / │
         /  │
lower──── 
```

$nonlinear_warning

Note: when composing linear systems with nonlinearities, it's often important to handle operating points correctly.
See [`ControlSystemsBase.offset`](@ref) for handling operating points.
"""
saturation(args...) = nonlinearity(Saturation(args...))
saturation(v::AbstractVector, args...) = nonlinearity(Saturation.(v, args...))
Base.show(io::IO, f::Saturation) = f.u == -f.l ? print(io, "saturation($(f.u))") : print(io, "saturation($(f.l), $(f.u))")

## Offset ======================================================================

struct Offset{T} <: Function
    o::T
end

(s::Offset)(x) = x .+ s.o

"""
    offset(val)

Create a constant-offset nonlinearity `x -> x + val`.

$nonlinear_warning

# Example:
To create a linear system that operates around operating point `y₀, u₀`, use
```julia
offset_sys = offset(y₀) * sys * offset(-u₀)
```
note the sign on the offset `u₀`. This ensures that `sys` operates in the coordinates
`Δu = u-u₀, Δy = y-y₀` and the inputs and outputs to the offset system are in their non-offset coordinate system. If the system is linearized around `x₀`, `y₀` is given by `C*x₀`. Additional information and an example is available here https://juliacontrol.github.io/ControlSystemsBase.jl/latest/lib/nonlinear/#Non-zero-operating-point
"""
offset(val::Number) = nonlinearity(Offset(val))
offset(v::AbstractVector) = nonlinearity(Offset.(v))

Base.show(io::IO, f::Offset) = print(io, "offset($(f.o))")



## Ratelimit ===================================================================

"""
    ratelimit(val; Tf)
    ratelimit(lower, upper; Tf)

Create a nonlinearity that limits the rate of change of a signal, roughly equivalent to `` 1/s ∘ sat ∘ s``. `Tf` controls the filter time constant on the derivative used to calculate the rate. 
$nonlinear_warning
"""
function ratelimit(args...; Tf)
    T = promote_type(typeof(Tf), typeof.(args)...)
    # Filtered derivative
    A = [0 1; -2/Tf^2 -2/Tf]
    B = [0; 1]
    C = [4/Tf^2 0]
    F = ss(A,B,C,0)

    integrator = tf([T(Tf), one(T)], [one(T), zero(T)])
    integrator*saturation(args...)*F
end



## DeadZone ====================================================================

struct DeadZone{T} <: Function
    l::T
    u::T
end
DeadZone(u) = DeadZone(-u, u)

function (s::DeadZone)(x)
    if x > s.u
        x - s.u
    elseif x < s.l
        x - s.l
    else
        zero(x)
    end
end

"""
    deadzone(val)
    deadzone(lower, upper)

Create a dead-zone nonlinearity.

```
       y▲
        │     /
        │    /
  lower │   /
─────|──┼──|───────► u
    /   │   upper
   /    │
  /     │
```

$nonlinear_warning

Note: when composing linear systems with nonlinearities, it's often important to handle operating points correctly.
See [`ControlSystemsBase.offset`](@ref) for handling operating points.
"""
deadzone(args...) = nonlinearity(DeadZone(args...))
deadzone(v::AbstractVector, args...) = nonlinearity(DeadZone.(v, args...))
Base.show(io::IO, f::DeadZone) = f.u == -f.l ? print(io, "deadzone($(f.u))") : print(io, "deadzone($(f.l), $(f.u))")