using ControlSystems
# NOTE: Glad, Ljung chap 9 contains several numerical examples that can be used as test cases

@testset "test_lqg" begin

function access_all_properties(G)
    pn = propertynames(G)
    for p in pn
        p === :stabilityrobustness && !ControlSystems.issiso(G.P) && continue
        @test_nowarn getproperty(G, p)
    end
end

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
access_all_properties(G)
gangoffourplot(G) # Test that it at least does not error
@test approxsetequal(eigvals(G.sysc.A), [ -31.6209+0.0im, -1.40629+0.0im, -15.9993+0.911174im, -15.9993-0.911174im, ], rtol = 1e-3)

qQ = 1
qR = 1
Q1 = 1000eye_(4)
Q2 = 1eye_(2)
R1 = 1eye_(6)
R2 = 1eye_(2)
N = eye_(6)
Gi = LQG(sys, Q1, Q2, R1, R2, qQ=qQ, qR=qR, integrator=true, ϵ=0.001, N=N)
access_all_properties(Gi)
gangoffourplot(Gi) # Test that it at least does not error
@test approxsetequal(eigvals(Gi.sysc.A), [-0.001, -0.001, -47.4832, -44.3442, -3.40255, -1.15355 ], rtol = 1e-2)

@test approxsetequal(eigvals(G.cl.A), [-1.0, -14.1774, -2.21811, -14.3206, -1.60615, -22.526, -1.0, -14.1774], rtol=1e-3)
@test approxsetequal(eigvals(G.T.A), [-22.526, -2.21811, -1.60615, -14.3206, -14.1774, -1.0, -1.0, -14.1774], rtol=1e-3)
@test approxsetequal(eigvals(G.S.A), [-22.526, -2.21811, -1.60615, -14.3206, -14.1774, -1.0, -1.0, -14.1774], rtol=1e-3)
@test approxsetequal(eigvals(G.CS.A), [-31.6209+0.0im, -1.40629+0.0im, -15.9993+0.911174im, -15.9993-0.911174im, -22.526+0.0im, -2.21811+0.0im, -1.60615+0.0im, -14.3206+0.0im, -14.1774+0.0im, -1.0+0.0im, -1.0+0.0im, -14.1774+0.0im], rtol=1e-3)
@test approxsetequal(eigvals(G.PS.A), [-1.0, -1.0, -3.0, -1.0, -22.526, -2.21811, -1.60615, -14.3206, -14.1774, -1.0, -1.0, -14.1774], rtol=1e-3)


@test approxsetequal(eigvals(Gi.cl.A), [-1.0, -44.7425, -44.8455, -2.23294, -4.28574, -2.06662, -0.109432, -1.31779, -0.78293, -1.0, -0.001, -0.001], rtol=1e-3)
@test approxsetequal(eigvals(Gi.T.A), [-44.7425, -44.8455, -4.28574, -0.109432, -2.23294, -2.06662, -1.31779, -0.78293, -1.0, -1.0], rtol=1e-3)
@test approxsetequal(eigvals(Gi.S.A), [-44.7425, -44.8455, -4.28574, -0.109432, -2.23294, -2.06662, -1.31779, -0.78293, -1.0, -1.0], rtol=1e-3)

@test eigvals(Gi.sysc.A) == eigvals(Gi.Fr.A)


## Test cases from Glad Ljung chap 9

# Aircrapt control

A = [-0.292 8.13 -201 9.77 0 -12.5 17.1
    -0.152  -2.54  0.561  -0.0004  0  107  7.68
    0.0364  -0.0678  -0.481  0.0012  0  4.67  -7.98
    0 1 0.0401 0 0 0 0
    0 0 1 0 0 0 0
    0 0 0 0 0 -20 0
    0 0 0 0 0 0 -20
]

B = [
    0 -2.15
    -31.7 0.0274
    0 1.48
    0 0
    0 0
    20 0
    0 20
]
# Example 9.1
M = [0 0 0 1 0 0 0
     0 0 0 0 1 0 0]
C = M

nx,ny,nu = size(A,1), size(C, 1), size(B,2)


Q1 = I(ny)
Q2 = I(nu)
R1 = I(nx)
R2 = I(ny)


sys = ss(A,B,C,0)
G = LQG(sys, Q1, Q2, R1, R2, M=M)
access_all_properties(G)
@test G.L ≈ [-0.0022 0.17 0.12 0.98 0.31 0.76 0.018
            -0.0038 0.028  -0.25 0.16  -0.95 0.072 0.10] rtol=0.01


Q2 = 0.1*I(nu)
G = LQG(sys, Q1, Q2, R1, R2, M=M)
@test G.L ≈ [-0.0036 0.39 0.21 3.11 0.63 1.54 0.046
            -0.0073 0.085 -0.78 0.57 -3.10 0.20 0.30] rtol=0.01

@test dcgain(G.cl) ≈ I



# Room temperature control
K1 = 1
K2 = 1
K3 = 1/8
K4 = 0.2
A = [
    -K3-K4 K4 0 0 0
    K1 -K1-K2 0 K2 0
    0 0 -0.01 -0.068 0
    0 0 1 0 0
    0 0 0 0 -0.0001
]
# @warn "I added the last -0.0001"

B = [
    K3
    0
    0
    0
    0
]

M = [0 1 0 0 0]
C = [
    0 1 0 0 1
    0 1 0 0 0
    0 0 0 1 0
]
N = [
    0 0
    0 0
    1 0
    0 0
    0 1
]
sys = ss(A,B,C,0)


Q1 = 1.0I(1)
Q2 = 0.01I(1)
R1 = diagm([10.0, 1])
# R1 = diagm([0, 0, 10.0, 1, 0])
R2 = diagm([0.0001, 0.1, 0.001])
G = LQG(sys, Q1, Q2, R1, R2, qQ = 0.000001, M=M, N=N)
access_all_properties(G)

# NOTE
# the kalman filter qutions will not work out until the noise model is incorporated

@test G.K ≈ [-0.0039  0.0003  0.0035
            0.0463  0.0027  0.9676
            8.5224  0.0148  99.7543
            1.1426  0.0097  14.1198
            99.9634  -0.0026 -0.8533] rtol=1e-2
@test G.L ≈ [2.8106 1.4071 6.8099 3.7874 0] rtol=1e-2

@test dcgain(G.cl)[] ≈ 1


## example 9.3

A = diagm([-1.0, -1, -2])
B = [
    1.0 0
    0.0 1
    0.0 1
]

C = [
    2.0 0 3
    1.0 1 0
]
M = C
sys = ss(A,B,C,0)


Q1 = I(2)
Q2 = 0.1I(2)
R1 = I(2)
R2 = 0.001I(2)
G = LQG(sys, Q1, Q2, R1, R2, integrator=true, M = M, measurement=true)
access_all_properties(G)

Lr = pinv(M * ((B * G.L[:, 1:3] - A) \ B))

@test Lr ≈ [
    4.24 -3.43
    0.12 4.82
] rtol = 0.01

# NOTE: Glad ljung does not have identity in the end here, but the eigenvalues of A-BL are identical for this L and their L
@test G.L ≈ [
    4.05 0.81 4.23 1 0 
    5.05 1.01 5.96 0 1
] rtol = 0.01

@test G.K ≈ [
    0 0
    0 0
    0 0
    31.62 0
    0 31.62
] rtol = 0.01

@test dcgain(G.cl) ≈ I

# Random system, compare with matlab output


A = [  0.154388  1.42524    0.914886  -0.313071
1.36189   0.30188   -0.767355   0.313058
-0.485766  0.926861  -0.755906  -0.309824
2.10783   0.146366  -1.65582    0.115643]

B = [
  0.625203  -1.87059
 -0.399363   0.516783
 -0.150984   0.161573
 -0.432394  -0.162178
]

C = [
  0.99517   -0.26967    0.255771  -0.0356366
 -0.358348   0.114812  -0.417828   0.983806
  0.45019   -0.166284   1.47825    2.21589
]

K = [ 1.22237  -0.114471  -1.10722
1.00941   0.445568  -1.91614]

P = ss(A,B,C,0)
C = ss(K)

SiA = [
 -0.3534     1.707    -3.451    -5.871
   1.674    0.1881    0.1495     1.257
 -0.3683    0.8871   -0.5026  -0.08282
   2.433   0.09358    -2.655    -1.636
]

SiB = [
  0.6252   -1.871
 -0.3994   0.5168
  -0.151   0.1616
 -0.4324  -0.1622
]

SiC = [
   -0.759    0.1587     1.276      2.61
  0.01776  -0.09757     2.761     3.844
]

SiD = [
    1      0
    0      1]

Si = ss(SiA, SiB, SiC, SiD)


@test isapprox(ControlSystems.input_sensitivity(P, C), Si, rtol=1e-3)


SoA = [
    -0.3534     1.707    -3.451    -5.871
      1.674    0.1881    0.1495     1.257
    -0.3683    0.8871   -0.5026  -0.08282
      2.433   0.09358    -2.655    -1.636
]

SoB = [
    1.124     0.905    -2.892
 -0.03348    -0.276     0.548
  0.02146  -0.08928    0.1424
   0.6922   0.02276   -0.7895
]

SoC = [
   0.9952   -0.2697    0.2558  -0.03564
  -0.3583    0.1148   -0.4178    0.9838
   0.4502   -0.1663     1.478     2.216 
]

SoD = I(3)

So = ss(SoA, SoB, SoC, SoD)


@test isapprox(ControlSystems.output_sensitivity(P, C), So, rtol=1e-3)



ATo  = [
      -0.3534     1.707    -3.451    -5.871
        1.674    0.1881    0.1495     1.257
      -0.3683    0.8871   -0.5026  -0.08282
        2.433   0.09358    -2.655    -1.636
]
 
BTo = [
    1.124     0.905    -2.892
    -0.03348    -0.276     0.548
    0.02146  -0.08928    0.1424
    0.6922   0.02276   -0.7895
    ]

CTo = [
    -0.9952   0.2697  -0.2558  0.03564
    0.3583  -0.1148   0.4178  -0.9838
    -0.4502   0.1663   -1.478   -2.216
]
To = ss(ATo, BTo, CTo, 0)

@test isapprox(ControlSystems.output_comp_sensitivity(P, C), To, rtol=1e-3)

ATi = [
      -0.3534     1.707    -3.451    -5.871
        1.674    0.1881    0.1495     1.257
      -0.3683    0.8871   -0.5026  -0.08282
        2.433   0.09358    -2.655    -1.636
]
 
BTi = [
    0.6252   -1.871
    -0.3994   0.5168
    -0.151   0.1616
    -0.4324  -0.1622
    ]

CTi = [
0.759   -0.1587    -1.276     -2.61
    -0.01776   0.09757    -2.761    -3.844
]
 
Ti = ss(ATi, BTi, CTi, 0)
 
@test isapprox(ControlSystems.input_comp_sensitivity(P, C), Ti, rtol=1e-3)


end