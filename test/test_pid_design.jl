module TestAnalysis
using CustomTest
using Base.Test
using ControlSystems


@test_approx_eq stabregionPID("exp(-sqrt(s))", logspace(-3,3), kd=1)[2][1] -1.022356911142034


P = tf(1,[1,1])^4
gangoffourplot(P,tf(1))
ωp = 0.8
kp,ki,C = loopshapingPI(P,ωp,phasemargin=60, doplot=true)
@test_approx_eq kp 0.82274734724854
@test_approx_eq ki -0.45472580052281375
