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



##

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
See [`ControlSystems.offset`](@ref) for handling operating points.
"""
deadzone(args...) = nonlinearity(DeadZone(args...))
deadzone(v::AbstractVector, args...) = nonlinearity(DeadZone.(v, args...))
Base.show(io::IO, f::DeadZone) = f.u == -f.l ? print(io, "deadzone($(f.u))") : print(io, "deadzone($(f.l), $(f.u))")