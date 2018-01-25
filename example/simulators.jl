abstract type AbstractSimulator end

# ============================================================================================

struct Simulator <: AbstractSimulator
    P::StateSpace
    f
    y
end
function Simulator(P::StateSpace,u = (t,x) -> 0)
    f = (t,x) -> P.A*x + P.B*u(t,x)
    y = (t,x) -> P.C*x
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

struct GainSchedulingSimulator <: AbstractSimulator
    s::Simulator
    f
    y
    r
    e
    controllers::Vector{StateSpace}
    conditions::Vector{T} where T <: Function
end
function GainSchedulingSimulator(P,ri,controllers::AbstractVector{Tu},conditions::AbstractVector{Tc}) where Tu <: LTISystem where Tc <: Function
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
        der[pinds] .+= s.P.B*u # Add input from active controller to system dynamics
        der
    end
    GainSchedulingSimulator(s,f,y,r,e,controllers,conditions)
end

# ============================================================================================

DiffEqBase.solve(s::AbstractSimulator, x0, tspan, args...; kwargs...) = solve(ODEProblem(s.f,x0,tspan),args...; kwargs...)
DiffEqBase.solve(s::GainSchedulingSimulator, x0, tspan, args...; kwargs...) = solve(ODEProblem(s.f,vcat(x0, zeros(sum(c->c.nx, s.controllers))),tspan),args...; kwargs...)

function jacobian(s::AbstractSimulator, x, u)
    ForwardDiff.jacobian
end
