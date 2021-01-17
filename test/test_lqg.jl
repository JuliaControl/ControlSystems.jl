
@testset "test_lqg" begin


w = exp10.(range(-5, stop=5, length=1000))
s = tf("s")
P = [1/(s+1) 2/(s+3); 1/(s+1) 1/(s+1)]
sys = ss(P)
sysmin = minreal(sys)
A,B,C,D = sys.A,sys.B,sys.C,sys.D
Am,Bm,Cm,Dm = sysmin.A,sysmin.B,sysmin.C,sysmin.D

@test approxsetequal(eigvals(Am), [-3,-1,-1])

Q1 = 100eye_(4)
Q2 = 1eye_(2)
R1 = 100eye_(4)
R2 = 1eye_(2)
G = LQG(sys, Q1, Q2, R1, R2)
gangoffourplot(G) # Test that it at least does not error
@test approxsetequal(eigvals(G.sysc.A), [ -31.6209+0.0im, -1.40629+0.0im, -15.9993+0.911174im, -15.9993-0.911174im, ], rtol = 1e-3)

qQ = 1
qR = 1
Q1 = 1000eye_(4)
Q2 = 1eye_(2)
R1 = 1eye_(6)
R2 = 1eye_(2)
Gi = LQG(sys, Q1, Q2, R1, R2, qQ=qQ, qR=qR, integrator=true)
gangoffourplot(Gi) # Test that it at least does not error
@test approxsetequal(eigvals(Gi.sysc.A), [0.0, 0.0, -47.4832, -44.3442, -3.40255, -1.15355 ], rtol = 1e-3)

@test approxsetequal(eigvals(G.cl.A), [-1.0, -14.1774, -2.21811, -14.3206, -1.60615, -22.526, -1.0, -14.1774], rtol=1e-3)
@test approxsetequal(eigvals(G.T.A), [-22.526, -2.21811, -1.60615, -14.3206, -14.1774, -1.0, -1.0, -14.1774], rtol=1e-3)
@test approxsetequal(eigvals(G.S.A), [-22.526, -2.21811, -1.60615, -14.3206, -14.1774, -1.0, -1.0, -14.1774], rtol=1e-3)
@test approxsetequal(eigvals(G.CS.A), [-31.6209+0.0im, -1.40629+0.0im, -15.9993+0.911174im, -15.9993-0.911174im, -22.526+0.0im, -2.21811+0.0im, -1.60615+0.0im, -14.3206+0.0im, -14.1774+0.0im, -1.0+0.0im, -1.0+0.0im, -14.1774+0.0im], rtol=1e-3)
@test approxsetequal(eigvals(G.PS.A), [-1.0, -1.0, -3.0, -1.0, -22.526, -2.21811, -1.60615, -14.3206, -14.1774, -1.0, -1.0, -14.1774], rtol=1e-3)


@test approxsetequal(eigvals(Gi.cl.A), [-1.0, -44.7425, -44.8455, -2.23294, -4.28574, -2.06662, -0.109432, -1.31779, -0.78293, -1.0, 0.0, 0.0], rtol=1e-3)
@test approxsetequal(eigvals(Gi.T.A), [-44.7425, -44.8455, -4.28574, -0.109432, -2.23294, -2.06662, -1.31779, -0.78293, -1.0, -1.0], rtol=1e-3)
@test approxsetequal(eigvals(Gi.S.A), [-44.7425, -44.8455, -4.28574, -0.109432, -2.23294, -2.06662, -1.31779, -0.78293, -1.0, -1.0], rtol=1e-3)



end
