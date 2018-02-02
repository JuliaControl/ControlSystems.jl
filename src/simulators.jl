abstract type AbstractSimulator end

# ============================================================================================

"""
    Simulator
# Fields:
    P::StateSpace
    f = (x,p,t) -> x
    y = (x,t)   -> y
"""
struct Simulator <: AbstractSimulator
    P::StateSpace
    f
    y
end


"""
    Simulator(P::StateSpace, u = (x,t) -> 0)

Used to simulate continuous-time systems. See function `?solve` for additional info.
# Usage:
```
using OrdinaryDiffEq
h              = 0.1
Tf             = 20
t              = 0:h:Tf
P              = ss(tf(1,[2,1])^2)
K              = 5
reference(x,t) = [1.]
s              = Simulator(P, reference)
os             = OutputFeedbackSimulator(P,reference,e->K*e)
x0             = [0.,0]
tspan          = (0.0,Tf)
sol            = solve(s, x0, tspan, Tsit5())
plot(t, s.y(sol, t)[:], lab="Open loop step response")
sol            = solve(os, x0, tspan, Tsit5())
plot!(t, os.y(sol, t)[:], lab="P controller K="*string(K))
```
"""
function Simulator(P::StateSpace, u = (x,t) -> 0)
    @assert iscontinuous(P) "Simulator only supports continuous-time system. See function `lsim` for simulation of discrete-time systems."
    f = (x,p,t) -> P.A*x + P.B*u(x,t)
    y(x,t) = P.C*x
    y(sol::ODESolution,t) = P.C*sol(t)
    Simulator(P, f, y)
end

# ============================================================================================

"""
    OutputFeedbackSimulator
# Fields:
    s = Simulator(P)
    e = (x,t) -> error
    u = (x,t) -> control signal
    s = Simulator(P,u)
    f = (x,p,t) -> x # Dynamics
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

Simulate system `P` with output error feedback `fb(r(x,t)-y)`

# Usage:
```
using OrdinaryDiffEq
h              = 0.1
Tf             = 20
t              = 0:h:Tf
P              = ss(tf(1,[2,1])^2)
K              = 5
reference(x,t) = [1.]
os             = OutputFeedbackSimulator(P,reference,e->K*e)
x0             = [0.,0]
tspan          = (0.0,Tf)
sol            = solve(os, x0, tspan, Tsit5())
plot(t, os.y(sol, t)', lab="P controller K="*string(K))
```
"""
function OutputFeedbackSimulator(P,r,fb)
    s = Simulator(P)
    e = (x,t) -> r(x,t) .- s.y(x,t)
    u = (x,t) -> fb(e(x,t))
    s = Simulator(P,u)
    f = (x,p,t) -> s.P.A*x + s.P.B*u(x,t)
    OutputFeedbackSimulator(s,f,s.y,r,u,e)
end

# ============================================================================================

"""
    StateFeedbackSimulator
# Fields:
    s = Simulator(P)
    r = (x,t) -> reference
    e = (x,t) -> r(x,t) .- x # Error
    u = (x,t) -> fb(e(x,t)) # Control signal
    s = Simulator(P,u)
    f = (x,p,t) -> s.P.A*x + s.P.B*u(x,t) # Dynamics
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
    e = (x,t) -> r(x,t) .- x
    u = (x,t) -> fb(e(x,t))
    s = Simulator(P,u)
    f = (x,p,t) -> s.P.A*x + s.P.B*u(x,t)
    StateFeedbackSimulator(s,f,s.y,r,u,e)
end

# ============================================================================================

"""
    GainSchedulingSimulator
# Fields:
    s::Simulator
    f = (x,p,t) -> x
    y = (x,t) -> y
    r = (x,t) -> reference
    e = (x,t) -> r(x,t) .- x # Error
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
    GainSchedulingSimulator(P,r,controllers::AbstractVector{LTISystem},conditions::AbstractVector{Function (x,y,r)->Bool}; inputfun=(u,t)->u)
`r` is a function (x,t)-> reference signal
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
reference(x,t) = 1.

@manipulate for kp1 = slider(linspace(0,20,50), label="kp1", value=5), ki1 = slider(linspace(0,20,50), label="ki1", value=1), kp2 = slider(linspace(0,20,50), label="kp2", value=5), ki2 = slider(linspace(0,20,50), label="ki2", value=1), th = linspace(0,1,10)
    s      = Simulator(Ps)
    os     = OutputFeedbackSimulator(Ps,reference,e->5e)
    x0     = [0.,0]
    tspan  = (0.0,Tf)
    sol    = solve(os, x0, tspan, Tsit5())
    plot(t, os.y(sol, t)[:], lab="P controller", show=false)


    controllers = [pid(kp=5, ki=1)]
    conditions  = [(x,y,r) -> y[1] < 0.5, (x,y,r) -> y[1] >= 0.5]
    tspan       = (0.0,Tf)
    x0          = [0.,0]
    gs          = GainSchedulingSimulator(Ps, reference, controllers, conditions)
    sol         = solve(gs, x0, tspan, Tsit5())
    plot!(t, gs.y(sol, t)[:], lab="PI controller", show=false)


    controllers = [pid(kp=kp1, ki=ki1), pid(kp=kp2, ki=ki2)]
    conditions  = [(x,y,r) -> y[1] < th, (x,y,r) -> y[1] >= th]
    tspan       = (0.0,Tf)
    x0          = [0.,0]
    gs          = GainSchedulingSimulator(Ps, reference, controllers, conditions)
    sol         = solve(gs, x0, tspan, Tsit5())
    plot!(t, gs.y(sol, t)[:], lab="Gain scheduled controller", ylims=(0,1.5), show=false)
    plot!([tspan...], [th, th], lab="Threshold", c=:black, l=:dash)
end
```
"""
function GainSchedulingSimulator(P,ri,controllers::AbstractVector{Tu},conditions::AbstractVector{Tc}; inputfun=(u,t)->u) where Tu <: StateSpace where Tc <: Function
    s = Simulator(P)
    pinds = 1:P.nx # Indices of plant-state derivative
    r = (x,t) -> ri(x[pinds],t)
    y(x,t) = s.y(x[pinds],t)
    y(sol::ODESolution,t) = P.C*sol(t)[pinds,:]
    e = (x,t) -> r(x,t) .- y(x,t)
    f = function(x,p,t)
        xyr = (x[pinds],y(x,t),r(x,t))
        index = findfirst(c->c(xyr...), conditions)
        @assert index > 0 "No condition returned true"
        der = similar(x)
        der[pinds] = P.A*x[pinds] # System dynamics
        et = e(x,t)
        ind = P.nx+1 # First index of currently updated controller
        for c in controllers
            c.nx == 0 && continue # Gain only, no states
            inds = ind:(ind+c.nx-1) # Controller indices
            der[inds] = c.A*x[inds] + c.B*et # Controller dynamics, driven by error
            ind += c.nx
        end
        c = controllers[index] # Active controller
        cind = P.nx + (index == 1 ? 1 : 1+sum(i->controllers[i].nx, 1:index-1)) # Get index of first controller state
        u = c.C*x[cind:(cind+c.nx-1)] + c.D*et # Form control signal
        der[pinds] .+= inputfun(P.B*u,t) # Add input from active controller to system dynamics
        der
    end
    GainSchedulingSimulator(s,f,y,r,e,controllers,conditions)
end

function GainSchedulingSimulator(P,ri,controllers::AbstractVector{Tu},conditions::AbstractVector{Tc}; inputfun=(u,t)->u) where Tu <: LTISystem where Tc <: Function
    GainSchedulingSimulator(P,ri,ss.(controllers),conditions; inputfun=inputfun)
end
# ============================================================================================

"""
    sol = solve(s::AbstractSimulator, x0, tspan, solver = Tsit5(), args...; kwargs...)
Simulate the system represented by `s` from initial state `x0` over time span `tspan = (t0,tf)`.
`args` and `kwargs` are sent to the `solve` function from `OrdinaryDiffEq`. The `solver` defaults to [`Tsit5()`](http://docs.juliadiffeq.org/stable/solvers/ode_solve.html)

See also `Simulator` `OutputFeedbackSimulator` `StateFeedbackSimulator` `GainSchedulingSimulator` `lsim`
"""
DiffEqBase.solve(s::AbstractSimulator, x0, tspan, solver = Tsit5(), args...; kwargs...) = solve(ODEProblem(s.f,x0,tspan), solver, args...; kwargs...)
DiffEqBase.solve(s::GainSchedulingSimulator, x0, tspan, solver = Tsit5(), args...; kwargs...) = solve(ODEProblem(s.f,vcat(x0, zeros(sum(c->c.nx, s.controllers))),tspan), solver, args...; kwargs...)
