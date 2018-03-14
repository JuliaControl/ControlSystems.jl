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
    @assert all(P.D .== 0) "Can not simulate systems with direct term D != 0"
    f = (dx,x,p,t) -> dx .= P.A*x .+ P.B*u(x,t)
    y(x,t) = P.C*x #.+ P.D*u(x,t)
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
    f = (dx,x,p,t) -> dx .= s.P.A*x .+ s.P.B*u(x,t)
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
    f = (dx,x,p,t) -> dx .= s.P.A*x .+ s.P.B*u(x,t)
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
    feedforward::Vector{Float64}
end
"""
    GainSchedulingSimulator(P,r,controllers::AbstractVector{LTISystem},conditions::AbstractVector{Function (x,y,r)->Bool}; inputfun=(u,t)->u, feedforward=zeros())
`r` is a function (x,t)-> reference signal
`controllers` is a list of LTIsystems representing the controllers
`conditions` is a list of function that accepts state, output and reference and returns a Bool indicating which controller is active. At least one function should return true for every input. If more than one condition is true, the index of the first condition in the list is chosen.
´feedforward´ is a vector of constants that are added to the control signal calculated by the active controller, i.e., `u = C*x + D*e + feedforward[active_index]`

For a usage example, see notebook in example folder.
```
"""
function GainSchedulingSimulator(P,ri,controllers::AbstractVector{Tu},conditions::AbstractVector{Tc};
    inputfun=(u,t)->u,
    feedforward=zeros(length(controllers))) where Tu <: StateSpace where Tc <: Function
    s = Simulator(P)
    pinds = 1:P.nx # Indices of plant-state derivative
    r = (x,t) -> ri(x[pinds],t)
    y(x,t) = s.y(x[pinds],t)
    y(sol::ODESolution,t) = P.C*sol(t)[pinds,:]
    e = (x,t) -> r(x,t) .- y(x,t)
    f = function(der,x,p,t)
        xyr = (x[pinds],y(x,t),r(x,t))
        index = findfirst(c->c(xyr...), conditions)
        @assert index > 0 "No condition returned true"
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
        der[pinds] = P.A*x[pinds] # System dynamics
        u = c.C*x[cind:(cind+c.nx-1)] .+ c.D*et .+ feedforward[index] # Form control signal
        size(inputfun(P.B*u,t)), size(der[pinds])
        der[pinds] .+= vec(inputfun(P.B*u,t)) # Add input from active controller to system dynamics
        der
    end
    GainSchedulingSimulator(s,f,y,r,e,controllers,conditions,feedforward)
end

function GainSchedulingSimulator(P,ri,controllers::AbstractVector{Tu},conditions::AbstractVector{Tc};
    inputfun=(u,t)->u,
    feedforward=zeros(length(controllers))) where Tu <: LTISystem where Tc <: Function
    GainSchedulingSimulator(P,ri,ss.(controllers),conditions; inputfun=inputfun, feedforward=feedforward)
end
# ============================================================================================

"""
    sol = solve(s::AbstractSimulator, x0, tspan,  args...; kwargs...)
Simulate the system represented by `s` from initial state `x0` over time span `tspan = (t0,tf)`.
`args` and `kwargs` are sent to the `solve` function from `OrdinaryDiffEq`, e.g., `solve(s, x0, tspan,  Tsit5(), reltol=1e-5)` solves the problem with solver [`Tsit5()`](http://docs.juliadiffeq.org/stable/solvers/ode_solve.html) and relative tolerance 1e-5.

See also `Simulator` `OutputFeedbackSimulator` `StateFeedbackSimulator` `GainSchedulingSimulator` `lsim`
"""
DiffEqBase.solve(s::AbstractSimulator, x0, tspan, solver=Tsit5(), args...; kwargs...) = solve(ODEProblem(s.f,x0,tspan), solver, args...; kwargs...)
DiffEqBase.solve(s::GainSchedulingSimulator, x0, tspan, solver=Tsit5(), args...; kwargs...) = solve(ODEProblem(s.f,vcat(x0, zeros(sum(c->c.nx, s.controllers))),tspan), solver, args...; kwargs...)
