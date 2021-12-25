@testset "test_analysis" begin
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
@test tzeros(ex_3) ≈ [0.3411639019140099 + 1.161541399997252im,
                             0.3411639019140099 - 1.161541399997252im,
                             0.9999999999999999 + 0.0im,
                             -0.6823278038280199 + 0.0im]
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
@test tzeros(ex_4) ≈ [-0.06467751189940692,-0.3680512036036696]

# Example 5
s = tf("s")
ex_5 = 1/s^15
@test tzeros(ex_5) == Float64[]
@test tzeros(ss(ex_5)) == Float64[]

# Example 6
A = [2 -1 0;
     0 0 0;
     -1 0 0]
B = [0; 0; 1]
C = [0 -1 0]
D = [0]
ex_6 = ss(A, B, C, D)
@test tzeros(ex_6) == Float64[]

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
# TODO: Improve the test for testing equalifity of sets of complex numbers
# i.e. simplify and handle doubles.
@test all(approxin(p,known_poles) for p in poles(ex_8)) && all(approxin2(p,poles(ex_8)) for p in known_poles)

ex_12 = ss(-3, 2, 1, 2)
@test poles(ex_12) ≈ [-3]

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
@test_broken true == false # order of poles changed below, should probably be consistent
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
@test markovparam(G, 0) == [0.0 0.0; 1.0 0.0]
@test markovparam(G, 1) == [0.0 3.0; -4.0 2.0]
@test markovparam(G, 2) == [0.0 -3.0; 20.0 -12.0]

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

# Example 5.5 from http://www.control.lth.se/media/Education/EngineeringProgram/FRTN10/2017/e05_both.pdf
G = [1/(s+2) -1/(s+2); 1/(s+2) (s+1)/(s+2)]
@test_broken length(poles(G)) == 1
@test length(tzeros(G)) == 1
@test_broken size(minreal(ss(G)).A) == (2,2)


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

end