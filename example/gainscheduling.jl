using ControlSystems, OrdinaryDiffEq, InteractNext
include("simulators.jl")

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
