@testset "test_pid_design" begin

# Test gof plot and loopshaping
P = tf(1,[1,1])^4
gangoffourplot(P,tf(1))
gangoffourplot(ss(P),ss(1))
ωp = 0.8
C, kp, ki = loopshapingPI(P, ωp, phasemargin=60, form=:parallel, doplot=true)
@test kp ≈ 0.82274734724854
@test ki ≈ -0.45472580052281375

# Test pid creation
# tf
@test pid(1.0, 1, 1) == tf(1) + tf(1,[1,0]) + tf([1,0],[1])
@test pid(1.0, Inf, 1) == tf(1) + tf([1, 0], [1])
# ss
@test tf(pid(1.0, 1, 0; state_space=true)) == tf(1) + tf(1,[1,0])

# Zero integral action
Tf = 0.01
@test tf(pid(2.0, 0; state_space=true, Tf)) ≈ minreal(pid(2.0, 0; state_space=false, Tf))

@test tf(pid(2.0, 0, 1; state_space=true, Tf)) ≈ minreal(pid(2.0, 0, 1; state_space=false, Tf))

# Test pidplots
C = pid(1.0, 1, 1) 
pidplots(C, :nyquist, :gof, :pz, :controller; params_p=[1.0, 2], params_i=[2, 3], grid=true) # Simply test that the functions runs and not errors
pidplots(C, :nyquist, :gof, :pz, :controller; params_p=1.0, params_i=[2, 3], grid=false)
leadlinkcurve()

# Test stabregionPID
stabregionPID(tf(1,[1,0]))
stabregionPID(s -> exp(-sqrt(s)), doplot=true)

w = exp10.(range(-3, stop=3, length=50))
kp, ki = stabregionPID(s->exp(-sqrt(s)), w; form=:parallel, kd=1)

for i in eachindex(w)
    w1 = w[i]
    @test exp(-sqrt(w1*im))*pid(kp[i],ki[i],1, form=:parallel)(w1*im)[] ≈ -1 atol=1e-2
end

# Test loopshapingPI
P = tf(1,[1, 1, 1])
C, kp, ki = loopshapingPI(P,10; phasemargin = 30, doplot = false)
_,_,_,pm = margin(P*C)
@test pm[] > 30

# Test placePI
P = tf(1,[1, 1])
C, Kp, Ti = placePI(P, 2, 0.7; form=:standard)
@test poles(feedback(P, C)) ≈ [-1.4 + √2.04im, -1.4 - √2.04im]
@test [Kp, Ti] ≈ [9/5, 9/20]

# Test internal functions convert_pidparams*
params = (2, 3, 0.5)
parallel_params = ControlSystems.convert_pidparams_from_standard(params..., :parallel)
@test parallel_params == (2, 2/3, 1)
@test ControlSystems.convert_pidparams_to_standard(parallel_params..., :parallel) == params
series_params = ControlSystems.convert_pidparams_from_standard(params..., :series)
@test series_params == ((3-sqrt(3))/3, (3-sqrt(3))/2, (3+sqrt(3))/2)
@test ControlSystems.convert_pidparams_to_standard(series_params..., :series) == params

# lead lag link
a = 1
M = 10
L = laglink(a, M)
@test freqresp(L, 1000a)[] ≈ 1+0im atol=1e-3
@test freqresp(L, 0)[] ≈ M+0im atol=1e-3

b = 1
K = 10
N = 5
L = leadlink(b, N, K)
@test freqresp(L, 10000a)[] ≈ K*N+0im atol=3e-2
@test freqresp(L, 0)[] ≈ K+0im atol=1e-3
@test abs(freqresp(L, b*√(N))[]) ≈ K*√(N) atol=1e-3

w = 1
L = leadlinkat(w, N, K)
@test abs(freqresp(L, w)[]) ≈ K*√(N) atol=1e-3


## Test loopshapingPID

Ms = 1.3
phasemargin = rad2deg(2asin(1/(2Ms)))
gm = Ms/(Ms-1)
P0 = tf(1, [1, 0, 0]) # loopshapingPI does not handle a double integrator
ωp = 1
form = :standard

# set rl = 1 to make the crossing point on the nyquistplot below easier to draw
C,kp,ki,kd = loopshapingPID(P, ωp; rl = 1, phasemargin)
@test kp >= 0
@test ki >= 0
@test kd >= 0

w = exp10.(LinRange(-0.1, 2, 200))
_,_,_,pm = margin(P*C, w)
@test pm[] > 0.99*phasemargin
@test rad2deg(angle(freqresp(P*C, ωp)[])) ≈ -180 + phasemargin atol=1e-2

# nyquistplot(P*C, w, unit_circle=true, Ms_circles=[Ms]); scatter!([cosd(-180+phasemargin)], [sind(-180+phasemargin)], lab="Specification point")

# Test one more system
P = let
    PA = [0.0 1.0; -1.2000000000000002e-6 -0.12000999999999999]
    PB = [0.0; 1.0;;]
    PC = [11.2 0.0]
    PD = [0.0;;]
    ss(PA, PB, PC, PD)
end
w = exp10.(LinRange(-2, 2, 400))
ωp = 4
C,kp,ki,kd = loopshapingPID(P, ωp; rl = 1, phasemargin)
@test kp >= 0
@test ki >= 0
@test kd >= 0

# nyquistplot(tf(P)*C, w, unit_circle=true, Ms_circles=[Ms]); scatter!([cosd(-180+phasemargin)], [sind(-180+phasemargin)], lab="Specification point")
end

