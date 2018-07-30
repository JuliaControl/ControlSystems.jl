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
x0             = [0.,0]
tspan          = (0.0,Tf)
sol            = solve(s, x0, tspan, Tsit5())
plot(t, s.y(sol, t)[:], lab="Open loop step response")
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


"""
    sol = solve(s::AbstractSimulator, x0, tspan,  args...; kwargs...)
Simulate the system represented by `s` from initial state `x0` over time span `tspan = (t0,tf)`.
`args` and `kwargs` are sent to the `solve` function from `OrdinaryDiffEq`, e.g., `solve(s, x0, tspan,  Tsit5(), reltol=1e-5)` solves the problem with solver [`Tsit5()`](http://docs.juliadiffeq.org/stable/solvers/ode_solve.html) and relative tolerance 1e-5.

See also `Simulator` `lsim`
"""
DiffEqBase.solve(s::AbstractSimulator, x0, tspan, solver=Tsit5(), args...; kwargs...) = solve(ODEProblem(s.f,x0,tspan), solver, args...; kwargs...)
