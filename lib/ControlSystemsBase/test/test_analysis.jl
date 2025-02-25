@test_throws MethodError poles(big(1.0)*ssrand(1,1,1)) # This errors before loading GenericSchur
using GenericSchur # Required to compute eigvals (in tzeros and poles) of a matrix with exotic element types
@testset "test_analysis" begin
es(x) = sort(x, by=LinearAlgebra.eigsortby)
## tzeros ##
# Examples from the Emami-Naeini & Van Dooren Paper
# Example 3
A = [0 1 0 0 0 0;
     0 0 1 0 0 0;
     0 0 0 0 0 0;
     0 0 0 0 1 0;
     0 0 0 0 0 1;
     0 0 0 0 0 0]
B = Matrix([0 0 1 0 0 0;
     0 0 0 0 0 1]')
C = [1 1 0 0 0 0;
     0 0 0 1 -1 0]
D = [1 0;
     1 0]

ex_3 = ss(A, B, C, D)
@test ControlSystemsBase.count_integrators(ex_3) == 6
@test es(tzeros(ex_3)) ≈ es([0.3411639019140099 + 1.161541399997252im,
                             0.3411639019140099 - 1.161541399997252im,
                             0.9999999999999999 + 0.0im,
                             -0.6823278038280199 + 0.0im])
# Example 4
A = [-0.129    0.0   0.396e-1  0.25e-1    0.191e-1;
     0.329e-2  0.0  -0.779e-4  0.122e-3  -0.621;
     0.718e-1  0.0  -0.1       0.887e-3  -0.385e1;
     0.411e-1  0.0   0.0      -0.822e-1   0.0;
     0.361e-3  0.0   0.35e-4   0.426e-4  -0.743e-1]
B = [0.0  0.139e-2;
     0.0  0.359e-4;
     0.0  -0.989e-2;
     0.249e-4  0.0;
     0.0  -0.534e-5]
C = [1 0 0 0 0;
     0 1 0 0 0]
D = zeros(2, 2)
ex_4 = ss(A, B, C, D)
@test es(tzeros(ex_4)) ≈ es([-0.06467751189940692,-0.3680512036036696])
@test ControlSystemsBase.count_integrators(ex_4) == 1

# Example 5
s = tf("s")
ex_5 = 1/s^15
@test tzeros(ex_5) == Float64[]
@test tzeros(ss(ex_5)) == Float64[]
@test ControlSystemsBase.count_integrators(ex_5) == 15

# Example 6
A = [2 -1 0;
     0 0 0;
     -1 0 0]
B = [0; 0; 1]
C = [0 -1 0]
D = [0]
ex_6 = ss(A, B, C, D)
@test tzeros(ex_6) == [2] # From paper: "Our algorithm will extract the singular part of S(A) and will yield a regular pencil containing the single zero at 2."
@test_broken tzeros(big(1.0)ex_6) == [2]
@test ControlSystemsBase.count_integrators(ex_6) == 2

@test ss(A, [0 0 1]', C, D) == ex_6

# Example 7
ex_7 = ss(zeros(2, 2), [0;1], [-1 0], [0])
@test tzeros(ex_7) == Float64[]

# Example 8
A = [-2 1 0 0 0 0;
     1 -2 1 0 1 -1;
     0 1 -2 1 0 0;
     0 0 1 -1 0 1;
     0 -1 0 0 0 0;
     0 1 0 -1 0 0]
B = [1; 0; 0; 0; 1; 0]
C = [0 0 0 1 0 0]
D = [0]
ex_8 = ss(A, B, C, D)
# TODO : there may be a way to improve the precision of this example.
@test tzeros(ex_8) ≈ [-1.0, -1.0] atol=1e-7
@test tzeros(big(1)ex_8) ≈ [-1.0, -1.0] atol=1e-7
@test ControlSystemsBase.count_integrators(ex_8) == 0

# Example 9
ex_9 = (s - 20)/s^15
@test tzeros(ex_9) ≈ [20.0]
@test tzeros(ss(ex_9)) ≈ [20.0]

# Example 11
A = [-2 -6 3 -7 6;
     0 -5 4 -4 8;
     0 2 0 2 -2;
     0 6 -3 5 -6;
     0 -2 2 -2 5]
B = Matrix([-2 -8 -3 1 -8;
     7 -5 0 5 0]')
C = [0 -1 2 -1 -1;
     1 1 1 0 -1;
     0 3 -2 3 -1]
D = [0 0;
     0 0;
     0 0]
ex_11 = ss(A, B, C, D)
@test tzeros(ex_11) ≈ [4.0, -3.0]
@test tzeros(big(1)ex_11) ≈ [4.0, -3.0]

# Test for multiple zeros, siso tf
s = tf("s")
sys = s*(s + 1)*(s^2 + 1)*(s - 3)/((s + 1)*(s + 4)*(s - 4))
@test tzeros(sys) ≈ [3.0, -1.0, im, -im, 0.0]

## POLE ##
@test poles(sys) ≈ [4.0, -4.0, -1.0]
@test poles([sys sys]) ≈ [4.0, -4.0, -1.0] # Issue #81
@test poles(ex_11) ≈ eigvals(ex_11.A)
@test poles([2/(s+1) 3/(s+2); 1/(s+1) 1/(s+1)]) ≈ [-1, -1, -2]


known_poles = [-3.383889568918823 + 0.000000000000000im
                            -2.199935841931115 + 0.000000000000000im
                            -0.624778101910111 + 1.343371895589931im
                            -0.624778101910111 - 1.343371895589931im
                            -0.083309192664918 + 0.487701968391972im
                            -0.083309192664918 - 0.487701968391972im]
approxin2(el,col) = any(el.≈col)
# Compares the computed poles with the expected poles
# TODO: Improve the test for testing equality of sets of complex numbers
# i.e. simplify and handle doubles.
@test all(approxin(p,known_poles) for p in poles(ex_8)) && all(approxin2(p,poles(ex_8)) for p in known_poles)

ex_12 = ss(-3, 2, 1, 2)
@test poles(ex_12) ≈ [-3]
@test ControlSystemsBase.count_integrators(ex_12) == 0

ex_13 = ss([-1 1; 0 -1], [0; 1], [1 0], 0)
@test poles(ex_13) ≈ [-1, -1]


## ZPKDATA ##

H = [tf(0) tf([3, 0],[1, 1, 10]) ; tf([1, 1],[1, 5]) tf([2],[1, 6])]
G = ss(H)
sol_z = vecarray(ComplexF64, 2, 2, ComplexF64[], ComplexF64[0.0 + 0.0im],
        ComplexF64[-1.0 + 0.0im], ComplexF64[])
sol_p = vecarray(ComplexF64, 2, 2, ComplexF64[], ComplexF64[-0.5 - 3.1224989991991996im,
        -0.5 + 3.1224989991991996im],
        ComplexF64[-5.0 + 0.0im], ComplexF64[-6.0 + 0.0im])
sol_k = [0.0 3.0; 1.0 2.0]
z, p, k = zpkdata(H)

@test_array_vecs_eps z sol_z 2*eps(Float64)
#@test_array_vecs_eps p sol_p 2*eps(Float64)
@test k == sol_k
z, p, k = zpkdata(G)
@test_array_vecs_eps z sol_z 10*eps(Float64)
@test_array_vecs_eps p sol_p 10*eps(Float64)
@test k == sol_k

## GAIN ## #Gain is confusing when referring to zpkdata. Test dcgain instead
@test [dcgain(H[1, 1]) dcgain(H[1, 2]); dcgain(H[2, 1]) dcgain(H[2, 2])] ≈ [0 0; 0.2 1/3]
@test [dcgain(G[1, 1]) dcgain(G[1, 2]); dcgain(G[2, 1]) dcgain(G[2, 2])] ≈ [0 0; 0.2 1/3]
@test dcgain(H[1, 1], 1e-6)[] == 0
@test dcgain(H[2, 1], 1e-6)[] ≈ 0.2 rtol=1e-5

## MARKOVPARAM ##
Hd = [tf(0, 1) tf([3, 0],[1, 1, 10], 1) ; tf([1, 1],[1, 5], 1) tf([2],[1, 6], 1)]
Gd = ss(Hd)
@test markovparam(Gd, 0) == [0.0 0.0; 1.0 0.0]
@test markovparam(Gd, 1) == [0.0 3.0; -4.0 2.0]
@test markovparam(Gd, 2) == [0.0 -3.0; 20.0 -12.0]

## DAMP ##
@test damp(sys)[1] ≈ [1.0, 4.0, 4.0]
@test damp(sys)[2] ≈ [1.0, -1.0, 1.0]

sysd = zpk([], [0.1+0.5im, 0.1-0.5im, 0.9], 1.0, 0.01)
@test damp(sysd)[1] ≈ [10.536051565782627, 152.96671271576292, 152.96671271576292]
@test damp(sysd)[2] ≈ [1.0, 0.4403159432698576, 0.4403159432698576]

damp_output = damp(ex_11)
@test damp_output[1] ≈ [1.0, 1.0, 2.0, 2.0, 3.0]
# THe order of the poles in ±1 and ±2 may come out in different order
@test damp_output[2][1:2] ≈ [-1.0, 1.0] || damp_output[2][1:2] ≈ [1.0, -1.0]
@test damp_output[2][3:4] ≈ [-1.0, 1.0] || damp_output[2][3:4] ≈ [1.0, -1.0]
@test damp_output[2][5] ≈ -1.0


## DAMPREPORT ##
s = tf("s")
sys = s*(s + 1)*(s^2 + 1)*(s - 3)/((s + 1)*(s + 4)*(s - 4))
@test sprint(dampreport, sys) == "|        Pole        |   Damping     |   Frequency   |   Frequency   | Time Constant |\n|                    |    Ratio      |   (rad/sec)   |     (Hz)      |     (sec)     |\n+--------------------+---------------+---------------+---------------+---------------+\n| -1                 |  1            |  1            |  0.159        |  1            |\n| +4                 |  -1           |  4            |  0.637        |  -0.25        |\n| -4                 |  1            |  4            |  0.637        |  0.25         |\n"

@test sprint(dampreport, ex_4) == "|        Pole        |   Damping     |   Frequency   |   Frequency   | Time Constant |\n|                    |    Ratio      |   (rad/sec)   |     (Hz)      |     (sec)     |\n+--------------------+---------------+---------------+---------------+---------------+\n| +0                 |  -1           |  0            |  0            |  -Inf         |\n| -0.0597 ± 0.0171im |  0.961        |  0.0621       |  0.00989      |  16.7         |\n| -0.0858            |  1            |  0.0858       |  0.0137       |  11.7         |\n| -0.18              |  1            |  0.18         |  0.0287       |  5.55         |\n"

@test sprint(dampreport, 1/(s+1+2im)/(s+2+3im)) == "|        Pole        |   Damping     |   Frequency   |   Frequency   | Time Constant |\n|                    |    Ratio      |   (rad/sec)   |     (Hz)      |     (sec)     |\n+--------------------+---------------+---------------+---------------+---------------+\n| -1            -2im |  0.447        |  2.24         |  0.356        |  1            |\n| -2            -3im |  0.555        |  3.61         |  0.574        |  0.5          |\n"

# Example 5.5 from https://www.control.lth.se/fileadmin/control/Education/EngineeringProgram/FRTN10/2019/e05_both.pdf
G = [1/(s+2) -1/(s+2); 1/(s+2) (s+1)/(s+2)]
@test_broken length(poles(G)) == 1 # The tf poles don't understand the cancellations
@test length(poles(ss(G, minimal=true))) == 1 # The ss version with minimal realization does
@test length(tzeros(G)) == 0 # tzeros converts to minimal ss relalization
@test minreal(ss(G)).A ≈ [-2]
@test ss(G, minimal=true).A ≈ [-2]


## MARGIN ##

w = exp10.(LinRange(-1, 2, 100))
P = tf(1,[1.0, 1])
ωgₘ, gₘ, ωϕₘ, ϕₘ = margin(P, w)
@test length.((ωgₘ, gₘ, ωϕₘ, ϕₘ)) == (1,1,1,1)
@test gₘ[] == Inf
@test ϕₘ[] == Inf


P = tf(1,[1.0, 1, 0])
ωgₘ, gₘ, ωϕₘ, ϕₘ = margin(P, w)
@test length.((ωgₘ, gₘ, ωϕₘ, ϕₘ)) == (1,1,1,1)
@test gₘ[] == Inf
@test ϕₘ[] ≥ 50
@test ωϕₘ[] ≈ 0.7871132039227572
marginplot(P, w)

ωgₘ, gₘ, ωϕₘ, ϕₘ = margin(P, w, allMargins=true)
@test length.((ωgₘ, gₘ, ωϕₘ, ϕₘ)) == (1,1,1,1)
@test gₘ[][] == Inf
@test ϕₘ[][] ≥ 50
@test ωϕₘ[][] ≈ 0.7871132039227572

## Delay margin
dm = delaymargin(P)[]
R = freqresp(P*delay(dm), ωϕₘ[][]-0.5:0.0001:ωϕₘ[][]+0.5)
@test minimum(abs, R .- (-1)) ≈ 0 atol=1e-3

P = 100DemoSystems.double_mass_model() * tf(1, [0.002, 1])
ωgₘ, gₘ, ωϕₘ, ϕₘ = margin(P, w, allMargins=true)
@test size(ϕₘ[1]) == (3,)
dm = delaymargin(P)[]
@test dm ≈ 0.03743336128009814 atol=1e-6
@test minimum(abs, R .- (-1)) ≈ 0 atol=1e-3

@test delaymargin(tf(0.1, [1, 1])) == Inf

# https://github.com/JuliaControl/ControlSystems.jl/issues/941
C = 0.6 + 30 * tf([1, 0], [1])
G = tf([0.04, 0.0001, 1.1], [1, 0.03, 254.9])
H = tf([0.25], [1, 1, 0.25])
L = C * G * H * delay(1)
m = margin(L)
Lw = freqresp(L, m[1][])[]
@test imag(Lw) ≈ 0 atol = 1e-6 # Test definition of gain margin
@test inv(-real(Lw)) ≈ m[2][] atol = 1e-6 # Test definition of gain margin

# https://github.com/JuliaControl/ControlSystems.jl/issues/961
P = tf(1,[5, 10.25, 6.25, 1])
w_180, gm, w_c, pm = margin(50P)
@test pm[] ≈ -35.1 rtol=1e-2

## Tricky case from https://www.reddit.com/r/ControlTheory/comments/1inhxsz/understanding_stability_in_highorder/
s = tf("s")
kpu = -10.593216768722073; kiu = -0.00063; t = 1000; tau = 180; a = 1/8.3738067325406132e-5;
kpd = 15.92190277847431; kid = 0.000790960718241793;
kpo = -10.39321676872207317; kio = -0.00063;
kpb = kpd; kib = kid;

C1 = (kpu + kiu/s)*(1/(t*s + 1))
C2 = (kpu + kiu/s)*(1/(t*s + 1))
C3 = (kpo + kio/s)*(1/(t*s + 1))
Cb = (kpb + kib/s)*(1/(t*s + 1))
OL = (ss(Cb)*ss(C1)*ss(C2)*ss(C3)*exp(-3*tau*s))/((C1 - a*s)*(C2 - a*s)*(C3 - a*s));

wgm, gm, ωϕₘ, ϕₘ = margin(OL; full=true, allMargins=true)
@test ϕₘ[][] ≈ -320 rtol=1e-2
for wgm in wgm[]
     @test mod(rad2deg(angle(freqresp(OL, wgm)[])), 360)-180 ≈ 0 atol=1e-1
end

nint = ControlSystemsBase.count_integrators(pade(OL, 2))
nintbig = ControlSystemsBase.count_integrators(big(1.0)*pade(OL, 2))
@test nint == nintbig # This one is very tricky and tests the threshold value of the eigval counting

@test ControlSystemsBase.count_integrators(pade(OL, 3)) == nintbig
@test ControlSystemsBase.count_integrators(pade(OL, 4)) == nintbig
@test ControlSystemsBase.count_integrators(pade(OL, 10)) == nintbig


# RGA
a = 10
P = ss([0 a; -a 0], I(2), [1 a; -a 1], 0)

w = exp10.(LinRange(-1, 2, 1000))
rgaplot(P, w)
rgaplot([P, 2P], w)

R = relative_gain_array(P, w)
@test maximum(abs, R) > 50 # Inf/NaN at w=10
@test minimum(abs, R) ≈ 1e-2 atol=1e-4

R = relative_gain_array(P, 9.99)
@test size(R) == size(P)
@test R[1,1] ≈ R[2,2] rtol=0.01
@test R[2,1] ≈ R[1,2] rtol=0.01
@test R[1,1] ≈ -R[1,2] rtol=0.01


# tests from https://arxiv.org/pdf/1805.10312.pdf
G = [
     4 2 2
     2 1 1
     2 1 1
]

R = relative_gain_array(G)
@test R ≈ fill(1/9, 3, 3)

A = [
     7 4 8
     7 2 5
     3 8 8
]
B = [
     21 16 16
     21 8 10
     9 32 16
]

R = [
     -2.47 -2.41 5.88
     3.29 0.94 -3.24
     0.18 2.47 -1.65
]

@test relative_gain_array(A) ≈ R atol=0.01
M = [A B]
@test relative_gain_array(M) ≈ 1/2 .* [R R] atol=0.01


## Gang of four
P = ssrand(2,3,2)
C = ssrand(3,2,2)

gof = gangoffour(P,C)
@test gof[1] == sensitivity(P,C)
@test gof[2] == G_PS(P,C)
@test gof[3] == G_CS(P,C)
@test gof[4] == comp_sensitivity(P,C)
@test_nowarn gangoffourplot(P, C)
@test_nowarn gangoffourplot([P, P], C)
@test_nowarn gangoffourplot(P, [C, C])

F = ssrand(2,2,2,proper=true)
@test_nowarn gangofseven(P, C, F)


## Approximate double integrator
P = let
     tempA = [
          0.0 0.0 1.0 0.0 0.0 0.0
          0.0 0.0 0.0 0.0 1.0 0.0
          -400.0 400.0 -0.4 -0.0 0.4 -0.0
          0.0 0.0 0.0 0.0 0.0 1.0
          10.0 -11.0 0.01 1.0 -0.011 0.001
          0.0 10.0 0.0 -10.0 0.01 -0.01
     ]
     tempB = [0.0 0.0; 0.0 0.0; 100.80166347371734 0.0; 0.0 0.0; 0.0 -0.001; 0.0 0.01]
     tempC = [0.0 0.0 0.0 1.0 0.0 0.0]
     tempD = [0.0 0.0]
     ss(tempA, tempB, tempC, tempD)
end
@test ControlSystemsBase.count_integrators(P) == 2

## Difficult test case for zeros
G = let
     tempA = [-0.6841991610512457 -0.0840213470263692 -0.0004269818661494616 -2.7625001165862086e-18; 2.081491248616774 0.0 0.0 8.404160870560225e-18; 0.0 24.837541148074962 0.12622006230897712 0.0; -1.2211265763794115e-14 -2.778983834717109e8 -1.4122312296634873e6 -4.930380657631326e-32]
     tempB = [-0.5316255605902501; 2.0811471051085637; -45.068824982602656; 5.042589978197361e8;;]
     tempC = [0.0 0.0 0.954929658551372 0.0]
     tempD = [0.0;;]
     ss(tempA, tempB, tempC, tempD)
end

@test length(tzeros(G)) == 3
@test es(tzeros(G)) ≈ es(tzeros(big(1)G))

end