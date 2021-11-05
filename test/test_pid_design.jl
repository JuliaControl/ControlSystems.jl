@testset "test_pid_design" begin
@test stabregionPID(s->exp(-sqrt(s)), exp10.(range(-3, stop=3, length=50)); kd=1)[2][1][1] ≈ -1.022356911142034

P = tf(1,[1,1])^4
gangoffourplot(P,tf(1))
gangoffourplot(ss(P),ss(1))
ωp = 0.8
C, kp, ki = loopshapingPI(P, ωp, phasemargin=60, form=:parallel, doplot=true)
@test kp ≈ 0.82274734724854
@test ki ≈ -0.45472580052281375

C = pid(1.0, 1, 1) 
@test C == tf(1) + tf(1,[1,0]) + tf([1,0],[1])
pidplots(C, :nyquist, :gof, :pz, :controller; params_p=[1.0, 2], params_i=[2, 3], grid=true) # Simply test that the functions runs and not errors
pidplots(C, :nyquist, :gof, :pz, :controller; params_p=1.0, params_i=[2, 3], grid=false)
leadlinkcurve()

stabregionPID(tf(1,[1,0]))
stabregionPID(s -> exp(-sqrt(s)))

P = tf(1,[1, 1, 1])
C, kp, ki = loopshapingPI(P,10; phasemargin = 30, doplot = false)
_,_,_,pm = margin(P*C)
@test pm[] > 30

P = tf(1,[1, 1])
C, kp, ki = placePI(P, 2, 0.7)
@test pole(feedback(P, C)) ≈ [-1.4 + √2.04im, -1.4 - √2.04im]
@test [kp, ki] ≈ [9/5, 9/20]

params = (2, 3, 0.5)
parallel_params = ControlSystems.convert_pidparams_from_standard(params..., :parallel)
@test parallel_params == (2, 2/3, 1)
@test ControlSystems.convert_pidparams_to_standard(parallel_params..., :parallel) == params
series_params = ControlSystems.convert_pidparams_from_standard(params..., :series)
@test series_params == ((3-sqrt(3))/3, (3-sqrt(3))/2, (3+sqrt(3))/2)
@test ControlSystems.convert_pidparams_to_standard(series_params..., :series) == params

end
