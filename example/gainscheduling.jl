using ControlSystems, OrdinaryDiffEq
# import Base.e
include("simulators.jl")
P             = tf(1,[2,1])^2
stepplot!([P, feedback(P,tf(5))],20)

h                     = 0.1
Tf                    = 20
t                     = 0:h:Tf
Ps                    = ss(P)
K                     = 5
reference(t,x)        = 1.
# y(t,x)                = Ps.C*x
# f(t,x)                = Ps.A*x + Ps.B*u(t,x)
#
# y(t,sol::ODESolution) = Ps.C*sol(t)
# e(t,x)                = reference(t) .- y(t,x)
# u(t,x)                = K .* e(t,x) # Must be an array



s       = Simulator(Ps)
os      = OutputeFeedbackSimulator(Ps,reference,e->5e)
x0      = [0.,0]
tspan   = (0.0,Tf)
sol     = solve(os, x0, tspan, Tsit5())

plot(t, os.y(t, sol)[:])
