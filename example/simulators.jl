abstract type AbstractSimulator end

# ============================================================================================

struct Simulator <: AbstractSimulator
    P::StateSpace
    f
    y
end
function Simulator(P::StateSpace,u = (t,x) -> 0)
    f = (t,x) -> P.A*x + P.B*u(t,x) # TODO: u can not appear here
    y = (t,x) -> Ps.C*x
    (T::typeof(y))(t,sol::ODESolution) = P.C*sol(t)
    Simulator(P, f, y)
end

# ============================================================================================

struct OutputFeedbackSimulator <: AbstractSimulator
    s::Simulator
    f
    y
    r
    u
    e
end
function OutputFeedbackSimulator(P,r,fb)
    s = Simulator(P)
    e = (t,x) -> r(t,x) .- s.y(t,x)
    u = (t,x) -> fb(e(t,x))
    s = Simulator(P,u)
    f = (t,x) -> s.P.A*x + s.P.B*u(t,x)
    OutputFeedbackSimulator(s,f,s.y,r,u,e)
end

# ============================================================================================

struct StateFeedbackSimulator <: AbstractSimulator
    s::Simulator
    f
    y
    r
    u
    e
end
function StateFeedbackSimulator(P,r,fb)
    s = Simulator(P)
    e = (t,x) -> r(t,x) .- x
    u = (t,x) -> fb(e(t,x))
    s = Simulator(P,u)
    f = (t,x) -> s.P.A*x + s.P.B*u(t,x)
    StateFeedbackSimulator(s,f,s.y,r,u,e)
end

# ============================================================================================

DiffEqBase.solve(s::AbstractSimulator, x0, tspan, args...; kwargs...) = solve(ODEProblem(s.f,x0,tspan),args...; kwargs...)
