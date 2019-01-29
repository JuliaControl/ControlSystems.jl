@testset "test_pid_design" begin
@test stabregionPID("exp(-sqrt(s))", exp10.(range(-3, stop=3, length=50)), kd=1)[2][1] ≈ -1.022356911142034

P = tf(1,[1,1])^4
gangoffourplot(P,tf(1))
ωp = 0.8
kp,ki,C = loopshapingPI(P,ωp,phasemargin=60, doplot=true)
@test kp ≈ 0.82274734724854
@test ki ≈ -0.45472580052281375
end
