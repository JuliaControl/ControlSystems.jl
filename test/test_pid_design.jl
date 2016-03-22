module TestAnalysis
using CustomTest
using ControlSystems


@test_approx_eq stabregionPID("exp(-sqrt(s))", logspace(-3,3), kd=1)[2][1] -1.022356911142034
