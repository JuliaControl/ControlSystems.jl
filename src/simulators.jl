abstract type AbstractSimulator end

# ============================================================================================

"""
    Simulator
# Fields:
    P::StateSpace
    f = (t,x) -> x
    y = (t,x) -> y
"""
struct Simulator <: AbstractSimulator
    P::StateSpace
    f
    y
end


"""
    Simulator(P::StateSpace, u = (t,x) -> 0)

Used to simulate continuous-time systems. See function `?solve` for additional info.
# Usage:
```
using OrdinaryDiffEq
h              = 0.1
Tf             = 20
t              = 0:h:Tf
P              = ss(tf(1,[2,1])^2)
K              = 5
reference(t,x) = [1.]
s              = Simulator(P, reference)
os             = OutputFeedbackSimulator(P,reference,e->K*e)
x0             = [0.,0]
tspan          = (0.0,Tf)
sol            = solve(s, x0, tspan, Tsit5())
plot(t, s.y(t, sol)[:], lab="Open loop step response")
sol            = solve(os, x0, tspan, Tsit5())
plot!(t, os.y(t, sol)[:], lab="P controller K="*string(K))
```
"""
function Simulator(P::StateSpace, u = (t,x) -> 0)
    @assert iscontinuous(P) "Simulator only supports continuous-time system. See function `lsim` for simulation of discrete-time systems."
    f = (t,x) -> P.A*x + P.B*u(t,x)
    y = (t,x) -> P.C*x
    (T::typeof(y))(t,sol::ODESolution) = P.C*sol(t)
    Simulator(P, f, y)
end

# ============================================================================================

"""
    OutputFeedbackSimulator
# Fields:
    s = Simulator(P)
    e = (t,x) -> error
    u = (t,x) -> control signal
    s = Simulator(P,u)
    f = (t,x) -> x # Dynamics
"""
struct OutputFeedbackSimulator <: AbstractSimulator
    s::Simulator
    f
    y
    r
    u
    e
end
"""
    OutputFeedbackSimulator(P,r,fb)

Simulate system `P` with output error feedback `fb(r(t,x)-y)`

# Usage:
```
using OrdinaryDiffEq
h              = 0.1
Tf             = 20
t              = 0:h:Tf
P              = ss(tf(1,[2,1])^2)
K              = 5
reference(t,x) = [1.]
os             = OutputFeedbackSimulator(P,reference,e->K*e)
x0             = [0.,0]
tspan          = (0.0,Tf)
sol            = solve(os, x0, tspan, Tsit5())
plot(t, os.y(t, sol)', lab="P controller K="*string(K))
```
"""
function OutputFeedbackSimulator(P,r,fb)
    s = Simulator(P)
    e = (t,x) -> r(t,x) .- s.y(t,x)
    u = (t,x) -> fb(e(t,x))
    s = Simulator(P,u)
    f = (t,x) -> s.P.A*x + s.P.B*u(t,x)
    OutputFeedbackSimulator(s,f,s.y,r,u,e)
end

# ============================================================================================

"""
    StateFeedbackSimulator
# Fields:
    s = Simulator(P)
    r = (t,x) -> reference
    e = (t,x) -> r(t,x) .- x # Error
    u = (t,x) -> fb(e(t,x)) # Control signal
    s = Simulator(P,u)
    f = (t,x) -> s.P.A*x + s.P.B*u(t,x) # Dynamics
"""
struct StateFeedbackSimulator <: AbstractSimulator
    s::Simulator
    f
    y
    r
    u
    e
end
"""
    StateFeedbackSimulator(P,r,fb)

"""
function StateFeedbackSimulator(P,r,fb)
    s = Simulator(P)
    e = (t,x) -> r(t,x) .- x
    u = (t,x) -> fb(e(t,x))
    s = Simulator(P,u)
    f = (t,x) -> s.P.A*x + s.P.B*u(t,x)
    StateFeedbackSimulator(s,f,s.y,r,u,e)
end

# ============================================================================================

"""
    GainSchedulingSimulator
# Fields:
    s::Simulator
    f = (t,x) -> x
    y = (t,x) -> y
    r = (t,x) -> reference
    e = (t,x) -> r(t,x) .- x # Error
    controllers::Vector{StateSpace}
    conditions::Vector{Function}
"""
struct GainSchedulingSimulator <: AbstractSimulator
    s::Simulator
    f
    y
    r
    e
    controllers::Vector{StateSpace}
    conditions::Vector{T} where T <: Function
end
"""
    GainSchedulingSimulator(P,r,controllers::AbstractVector{LTISystem},conditions::AbstractVector{Function (x,y,r)->Bool})
`r` is a function (t,x)-> reference signal
`controllers` is a list of LTIsystems representing the controllers
`conditions` is a list of function that accepts state, output and reference and returns a Bool indicating which controller is active. At least one function should return true for every input. If more than one condition is true, the index of the first condition in the list is chosen.

# Example:
This example relies on `InteractNext.jl` for interactivity. If InteracNext is not installed and configured, manually select values for `kp1,ki1,kp2,ki2` and remove the `@manipulate` block.
```julia
using ControlSystems, OrdinaryDiffEq, InteractNext

P              = tf(1,[2,1])^2
h              = 0.1
Tf             = 20
t              = 0:h:Tf
Ps             = ss(P)
K              = 5
reference(t,x) = 1.

@manipulate for kp1 = slider(linspace(0,20,50), label="kp1", value=5), ki1 = slider(linspace(0,20,50), label="ki1", value=1), kp2 = slider(linspace(0,20,50), label="kp2", value=5), ki2 = slider(linspace(0,20,50), label="ki2", value=1), th = linspace(0,1,10)
    s      = Simulator(Ps)
    os     = OutputFeedbackSimulator(Ps,reference,e->5e)
    x0     = [0.,0]
    tspan  = (0.0,Tf)
    sol    = solve(os, x0, tspan, Tsit5())
    plot(t, os.y(t, sol)[:], lab="P controller", show=false)


    controllers = [pid(kp=5, ki=1)]
    conditions  = [(x,y,r) -> y[1] < 0.5, (x,y,r) -> y[1] >= 0.5]
    tspan       = (0.0,Tf)
    x0          = [0.,0]
    gs          = GainSchedulingSimulator(Ps, reference, controllers, conditions)
    sol         = solve(gs, x0, tspan, Tsit5())
    plot!(t, gs.y(t, sol)[:], lab="PI controller", show=false)


    controllers = [pid(kp=kp1, ki=ki1), pid(kp=kp2, ki=ki2)]
    conditions  = [(x,y,r) -> y[1] < th, (x,y,r) -> y[1] >= th]
    tspan       = (0.0,Tf)
    x0          = [0.,0]
    gs          = GainSchedulingSimulator(Ps, reference, controllers, conditions)
    sol         = solve(gs, x0, tspan, Tsit5())
    plot!(t, gs.y(t, sol)[:], lab="Gain scheduled controller", ylims=(0,1.5), show=false)
    plot!([tspan...], [th, th], lab="Threshold", c=:black, l=:dash)
end
```
"""
function GainSchedulingSimulator(P,ri,controllers::AbstractVector{Tu},conditions::AbstractVector{Tc}; inputfun=identity) where Tu <: LTISystem where Tc <: Function
    controllers = ss.(controllers)
    s = Simulator(P)
    pinds = 1:P.nx # Indices of plant-state derivative
    r = (t,x) -> ri(t,x[pinds])
    y = (t,x) -> s.y(t,x[pinds])
    (T::typeof(y))(t,sol::ODESolution) = P.C*sol(t)[pinds,:]
    e = (t,x) -> r(t,x) .- y(t,x)
    f = function(t,x)
        xyr = (x[pinds],y(t,x),r(t,x))
        index = length(controllers) > 1 ? findfirst(c->c(xyr...), conditions) : 1
        der = similar(x)
        der[pinds] = P.A*x[pinds] # System dynamics
        et = e(t,x)
        ind = P.nx+1 # First index of currently updated controller
        for c in controllers
            c.nx == 0 && continue
            inds = ind:ind+c.nx-1 # Controller indices
            der[inds] = c.A*x[inds] + c.B*et # Controller state derivatives
            ind += c.nx
        end
        c = controllers[index] # Active controller
        cind = P.nx + (index == 1 ? 1 : sum(c->c.nx, controllers[1:index-1])) # Get index of first controller state
        u = c.C*x[cind:cind+c.nx-1] + c.D*et # Form control signal
        der[pinds] .+= inputfun(s.P.B*u) # Add input from active controller to system dynamics
        der
    end
    GainSchedulingSimulator(s,f,y,r,e,controllers,conditions)
end

# ============================================================================================

"""
    sol = solve(s::AbstractSimulator, x0, tspan, args...; kwargs...)
Simulate the system represented by `s` from initial state `x0` over time span `tspan = (t0,tf)`.
`args` and `kwargs` are sent to the `solve` function from `OrdinaryDiffEq`

See also `Simulator` `OutputFeedbackSimulator` `StateFeedbackSimulator` `GainSchedulingSimulator` `lsim`
"""
DiffEqBase.solve(s::AbstractSimulator, x0, tspan, args...; kwargs...) = solve(ODEProblem(s.f,x0,tspan),args...; kwargs...)
DiffEqBase.solve(s::GainSchedulingSimulator, x0, tspan, args...; kwargs...) = solve(ODEProblem(s.f,vcat(x0, zeros(sum(c->c.nx, s.controllers))),tspan),args...; kwargs...)

function jacobian(s::AbstractSimulator, x, u)
    error("Not yet implemented")
end
