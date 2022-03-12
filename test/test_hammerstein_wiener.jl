using ControlSystems
using ControlSystems: HammersteinWienerSystem

# @testset "test_nonlinearity_system" begin

# For simplicity, equality of HammersteinWienerSystems are tested over a finite set of frequencies
ω = 0.0:8

@test typeof(promote(nonlinearity(abs2), ss(1))[1]) == HammersteinWienerSystem{Float64}

@test sprint(show, ss(1,1,1,1)*nonlinearity(abs2)) == "HammersteinWienerSystem{Float64}\n\nP: StateSpace{Continuous, Float64}\nA = \n 1.0\nB = \n 0.0  1.0\nC = \n 1.0\n 0.0\nD = \n 0.0  1.0\n 1.0  0.0\n\nContinuous-time state-space model\n\nNonlinearities: Any[abs2]"


P1 = HammersteinWienerSystem(ss(-1.0, 1, 1, 0))
P2 = HammersteinWienerSystem(ss(-2.0, -1, 1, 1)) # (s+1)/(s+2)

# Equality
@test P1 == deepcopy(P1)
@test P1 != deepcopy(P2)


# Random conversions
s = tf("s")
sys1 = HammersteinWienerSystem(1.0/s)
@test sys1.P.A == sys1.P.D == fill(0,1,1)
@test sys1.P.B*sys1.P.C == fill(1,1,1)
@test sys1.Tau == []

## Multiple nonlinearitys
G = ss(1.0) + nonlinearity(abs2) + nonlinearity(abs2)

s11 = feedback(ss(1/s), nonlinearity(abs2))
s12 = ss(1/s)
s21 = HammersteinWienerSystem(ss(1/(s+1)))
s22 = ss(1/(s+10))

s1 = [s11 s12]
s2 = [s21 s22]

f1 = [s1;s2]
f2 = [s11 s12;
      s21 s22]


@test propertynames(nonlinearity(abs2)) == (:P, :Tau, :nu, :ny)

# Test step
println("Simulating first nonlinearity system:")
@test step(nonlinearity(identity)*tf(1,[1.,1]), 5).y ≈ step(tf(1,[1.,1]), 5).y rtol=1e-3
@test step(nonlinearity(abs2)*tf(1,[1.,1]), 5).y ≈ abs2.(step(tf(1,[1.,1]), 5).y) rtol=1e-3
@test step(tf(1,[1.,1])*nonlinearity(x->0.5x), 5).y ≈ 0.5step(tf(1,[1.,1]), 5).y rtol=1e-3

@time step(nonlinearity(abs2)*tf(1,[1.,1]))
@time step(nonlinearity(abs2)*tf(1,[1,1]))
@test step(nonlinearity(abs2)*tf(1,[1,1])) isa ControlSystems.SimResult

res = step([s11;s12], 10)
@test size(res.u,1) == 1
@time y1, t1, x1 = res
@time @test y1[2:2,:] ≈ step(s12, t1)[1] rtol = 1e-14

res = step([s11 s12], 10)
@test size(res.u,1) == 2
@test size(res.u,3) == 2

t = 0.0:0.1:10

y1, t1, x1 = step(HammersteinWienerSystem([1.0/s 2/s; 3/s 4/s]), t)
y2, t2, x2 = step([1.0/s 2/s; 3/s 4/s], t)
@test y1 ≈ y2 rtol=1e-14
@test size(x1,2) == length(t)
@test size(x1,3) == 2


## Test nonlinearity with D22 term
t = 0:0.01:4

sys = nonlinearity(abs2)

y, t, x = step(sys, t)
@test y[:] ≈ ones(length(t)) atol = 1e-12
@test size(x) == (0,401)


## Test that nonlinearities can be composed 
sys = nonlinearity(abs2)*nonlinearity(abs2)*1/s
res = step(sys, t)

sys2 = nonlinearity(abs2 ∘ abs2)*1/s
res2 = step(sys2, t)

@test res.y ≈ res2.y rtol=1e-3

@error "jag kom inte längre hän hit. Är problem med att stoppa två olinjäriteter i rad, jag stoppade in assrt D22=0 för att undvika algebraisk loop i solvern, men det innebär att man måste fånga fallet ovan och göra kompositionen av olinjäriteterna"

y_sol = [zeros(200);0:0.01:2]'

@test maximum(abs,y-y_sol) < 1e-13
@test maximum(abs,x-collect(0:0.01:4)') < 1e-15

# TODO For some reason really bad accuracy here
# Looks like a lag in time
sys = 1/s*nonlinearity(abs2)*nonlinearity(abs2)

y, t, x = step(sys, t)
@test maximum(abs,y-y_sol) < 1e-5
@test maximum(abs,x-y_sol) < 1e-5

t = 0:0.001:0.1
y, t, x = step(sys, t)
@test length(y) == length(t)


##
G = tf(1.0, [1,1])
nl = nonlinearity(x->clamp(x, -0.7, 0.7))
Gnl = G*nl
C = tf(1)
Cnl = nl*C
L = G*C
Lnl = G*Cnl

plot(step([feedback(L); feedback(C,G)], 5), lab=["Linear y" "Linear u"])
plot!(step([feedback(Lnl); feedback(Cnl,G)], 5), lab=["Nonlinear y" "Nonlinear u"])


# end
