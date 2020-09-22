@testset "test_timeevol" begin


@test ControlSystems.common_timeevol(Discrete(1.5), Discrete(1.5)) == Discrete(1.5)
@test ControlSystems.common_timeevol(Continuous(), Continuous()) == Continuous()

@test_throws ErrorException ControlSystems.common_timeevol(Continuous(), Discrete(1.5))
@test_throws ErrorException ControlSystems.common_timeevol(Discrete(1.0), Discrete(1.5))


sys = ss(-1, 1, 1, 0)
iscontinuous(sys) == true
isdiscrete(sys) == false

sysd = ss(-1, 1, 1, 0, 0.1)
iscontinuous(sysd) == false
isdiscrete(sysd) == true


end
