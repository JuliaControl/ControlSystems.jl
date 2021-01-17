@testset "test_pid_design" begin
@test stabregionPID(s->exp(-sqrt(s)), exp10.(range(-3, stop=3, length=50)), kd=1)[2][1] ≈ -1.022356911142034

P = tf(1,[1,1])^4
gangoffourplot(P,tf(1))
gangoffourplot(ss(P),ss(1))
ωp = 0.8
kp,ki,C = loopshapingPI(P,ωp,phasemargin=60, doplot=true)
@test kp ≈ 0.82274734724854
@test ki ≈ -0.45472580052281375

C = pid(;kp = 1, ki = 1, kd=1) 
@test C == tf(1) + tf(1,[1,0]) + tf([1,0],[1])
pidplots(C, :nyquist,:gof,:pz,:controller, grid=true) # Simply test that the functions runs and not errors
pidplots(C, :nyquist,:gof,:pz,:controller, grid=false)
leadlinkcurve()

stabregionPID(tf(1,[1,0]))
stabregionPID(s -> exp(-sqrt(s)))

P = tf(1,[1, 1, 1])
kp,ki,C = loopshapingPI(P,10; phasemargin = 30, doplot = false)
_,_,_,pm = margin(P*C)
@test pm[] > 30

end
