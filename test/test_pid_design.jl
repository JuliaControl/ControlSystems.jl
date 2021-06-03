@testset "test_pid_design" begin
@test stabregionPID(s->exp(-sqrt(s)), exp10.(range(-3, stop=3, length=50)); kd=1, form=:parallel)[2][:Kp][1] ≈ -1.022356911142034

P = tf(1,[1,1])^4
gangoffourplot(P,tf(1))
gangoffourplot(ss(P),ss(1))
ωp = 0.8
C, params = loopshapingPI(P, ωp, phasemargin=60, form=:parallel, doplot=true)
@test params.Kp ≈ 0.82274734724854
@test params.Ki ≈ -0.45472580052281375

C = pid(TransferFunction; Kp=1.0, Ki=1, Kd=1) 
@test C == tf(1) + tf(1,[1,0]) + tf([1,0],[1])
pidplots(C, :nyquist, :gof, :pz, :controller; Kp=[1.0, 2], Ki=[2, 3], grid=true) # Simply test that the functions runs and not errors
pidplots(C, :nyquist, :gof, :pz, :controller; Kp=1.0, Ki=[2, 3], grid=false)
leadlinkcurve()

stabregionPID(tf(1,[1,0]))
stabregionPID(s -> exp(-sqrt(s)))

P = tf(1,[1, 1, 1])
C, params = loopshapingPI(P,10; phasemargin = 30, doplot = false)
_,_,_,pm = margin(P*C)
@test pm[] > 30

P = tf(1,[1, 1])
C, params = placePI(P, 2, 0.7)
@test pole(feedback(P, C)) ≈ [-1.4 + √2.04im, -1.4 - √2.04im]
@test [params.Kp, params.Ti] ≈ [9/5, 9/20]

@test all(values(ControlSystems.convert_pid_params(:standard; Kp=2.0, Ki=3, Kd=4)) .≈ (2, 2/3, 2))
@test all(values(ControlSystems.convert_pid_params(:series; Kp=2.0, Ti=4, Td=1)) .≈ (1, 2, 2))
@test all(values(ControlSystems.convert_pid_params(:parallel; Kc=2.0, τi=3, τd=4)) .≈ (14/3, 2/3, 8))

end
