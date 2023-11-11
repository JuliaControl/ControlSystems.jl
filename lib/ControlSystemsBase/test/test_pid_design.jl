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
@test pid(1.0, 0, 1) == tf(1) + tf([1, 0], [1])
@test pid(0.0, 1, 1; form=:parallel) == tf(0) + tf(1,[1,0]) + tf([1,0],[1])
# ss
@test tf(pid(1.0, 1, 0; state_space=true)) == tf(1) + tf(1,[1,0])

# Discrete
@test_throws ArgumentError pid(1.0, 1, 1, Ts=0.1)
@test_nowarn pid(1.0, 1, 0, Ts=0.1) # Without D term it's okay

Cd = pid(1.0, 1, 1, Ts=0.1, Tf=0.1)
@test isdiscrete(Cd)
@test Cd.Ts == 0.1

Cd2 = pid(1.0, 1, 1, Ts=0.01, Tf=0.1, state_space=true)
@test isdiscrete(Cd2)
@test Cd2 isa StateSpace
@test Cd2.Ts == 0.01
@test tzeros(Cd2) != tzeros(Cd) # Different sample rates

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

C, kp, ki, fig, CF = loopshapingPI(P,1; phasemargin = 45, doplot = true, Tf = 0.1)
_,_,_,pm = margin(P*CF)
@test pm[] > 45*0.99


# Test placePI
P = tf(1,[1, 1])
C, Kp, Ti = placePI(P, 2, 0.7; form=:standard)
@test poles(feedback(P, C)) ≈ [-1.4 + √2.04im, -1.4 - √2.04im]
@test [Kp, Ti] ≈ [9/5, 9/20]

# Test internal functions convert_pidparams*
params = (4, 3, 0.5)
standard_params = ControlSystemsBase.convert_pidparams_from_parallel(params..., :standard)
@test standard_params == (4, 4/3, 0.5/4)
@test ControlSystemsBase.convert_pidparams_to_parallel(standard_params..., :standard) == params
series_params = ControlSystemsBase.convert_pidparams_from_parallel(params..., :series)
@test series_params == ((4-sqrt(10))/2, (4-sqrt(10))/6, (4+sqrt(10))/6)
@test all(ControlSystemsBase.convert_pidparams_to_parallel(series_params..., :series) .≈ params)

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

Mt = Ms = 1.3
phasemargin = rad2deg(2asin(1/(2Ms)))
gm = Ms/(Ms-1)
P = tf(1, [1, 0, 0]) # loopshapingPI does not handle a double integrator
ω = 1
form = :standard

# set rl = 1 to make the crossing point on the nyquistplot below easier to draw
Mt = 1.3
ϕt = 75
C,kp,ki,kd = loopshapingPID(P, ω; Mt, ϕt)
@test kp >= 0
@test ki >= 0
@test kd >= 0

T = comp_sensitivity(P, C)
mt, wt = hinfnorm(T)
@test mt ≈ Mt rtol=1e-2
@test wt ≈ ω  rtol=1e-2




# Test one more system
P = let
    PA = [0.0 1.0; -1.2000000000000002e-6 -0.12000999999999999]
    PB = [0.0; 1.0;;]
    PC = [11.2 0.0]
    PD = [0.0;;]
    ss(PA, PB, PC, PD)
end

ω = 4
C,kp,ki,kd,fig = loopshapingPID(tf(P), ω; Mt, ϕt, doplot=true)
@test kp >= 0
@test ki >= 0
@test kd >= 0

T = comp_sensitivity(tf(P), C)
mt, wt = hinfnorm(T)
@test mt ≈ Mt rtol=1e-2
@test wt ≈ ω  rtol=1e-2


##
P = ControlSystemsBase.DemoSystems.resonant()
ω = 5
Mt = 1.3
ϕt = 75
C,kp,ki,kd,fig = loopshapingPID(tf(P), ω; Mt, ϕt, doplot=true)
@test kp >= 0
@test ki >= 0
@test kd >= 0

T = comp_sensitivity(tf(P), C)
mt, wt = hinfnorm(T)
@test mt ≈ Mt rtol=1e-2
@test wt ≈ ω  rtol=1e-2

##
P = ControlSystemsBase.DemoSystems.double_mass_model()
ω = 100
for ω = [1, 100]
    Mt = 1.3
    ϕt = 75
    C,kp,ki,kd,fig = loopshapingPID(tf(P), ω; Mt, ϕt, doplot=true)
    @test kp >= 0
    @test ki >= 0
    @test kd >= 0

    T = comp_sensitivity(tf(P), C)
    mt, wt = hinfnorm(T)
    @test mt ≈ Mt rtol=1e-2
    @test wt ≈ ω  rtol=1e-2
end

##
P = ControlSystemsBase.DemoSystems.sotd()
ω = 0.5
Mt = 1.1
ϕt = 10
C,kp,ki,kd,fig,CF = loopshapingPID((P), ω; Mt, ϕt, doplot=true)
@test kp >= 0
@test ki >= 0
@test kd >= 0

T = comp_sensitivity((P), CF)
w = exp10.(LinRange(-2, 2, 2000))
mag = ControlSystemsBase.bode(T, w)[1]
mt = maximum(mag)
wt = w[findmax(mag[:])[2]]
@test mt ≈ Mt rtol=1e-2
@test wt ≈ ω  rtol=1e-2

end

