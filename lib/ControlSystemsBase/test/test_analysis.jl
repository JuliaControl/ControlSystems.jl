@test_throws MethodError poles(big(1.0)*ssrand(1,1,1)) # This errors before loading GenericSchur
using GenericSchur # Required to compute eigvals (in tzeros and poles) of a matrix with exotic element types

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
                             -0.6823278038280199 + 0.0im]) atol=1e-6
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
@test tzeros(big(1.0)ex_6) == [2]
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
@test tzeros(ex_8) ≈ [-1.0, -1.0] atol=1e-6
@test tzeros(big(1)ex_8) ≈ [-1.0, -1.0] atol=1e-12
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
@test tzeros(ex_11) ≈ [4.0, -3.0] atol=1e-5
@test tzeros(big(1)ex_11) ≈ [4.0, -3.0]

# Test for multiple zeros, siso tf
s = tf("s")
sys = s*(s + 1)*(s^2 + 1)*(s - 3)/((s + 1)*(s + 4)*(s - 4))
@test tzeros(sys) ≈ [3.0, -1.0, im, -im, 0.0]

## POLE ##
@test es(poles(sys)) ≈ es([4.0, -4.0, -1.0])
@test es(poles(sys, atol=1e-6)) ≈ es([4.0, -4.0]) # This cancels the -1 pole/zero
@test es(poles([sys sys], atol=1e-6)) ≈ es([4.0, -4.0]) # Issue #81
@test poles(zpk([sys sys])) ≈ [4.0, -4.0, -1.0] # This does not cancel the -1 pole/zero
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
@test length(poles(G)) == 1 # The tf poles do understand the cancellations
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

##
G = let # From https://juliacomputing.github.io/JuliaSimControl.jl/v0.10/examples/mtk_disturbance_modeling/#Controller-reduction
     marginplottestA = [-0.8778418152658126 -0.5578161193611026 0.6099186362718944 -0.739076573488832 0.0 -0.7985033210421671 -0.6437387074246465 -1.09990118296936 -0.9415890954623862 -0.6105710677920562 -0.40609522036480955 -0.02167670496721598 -0.0922028349547116 -0.3659595949709856 0.0 0.0 0.0 0.0; -0.5578161193611026 -0.8220564196754826 -0.4879223680823921 0.8803386749291978 0.0 -0.5305709989906725 -0.42546859154084343 -0.7310995753661961 -0.6252826290802145 -0.4179592857573079 -0.274638775913742 -0.00903355930082124 -0.05630074104219519 -0.24851140114082376 0.0 0.0 0.0 0.0; -12.772894974281625 7.422755287471556 -8.893245512259899 -0.5983738877209197 4.4676642940721933e-7 -2.901639071283097 -2.3556062284249624 -3.998192111860942 -3.428202120843976 -2.1327661674537324 -1.4414918123884637 -0.11684112979233667 -0.37069980312317624 -1.2936533657321307 0.0 0.0 0.0 0.0; 9.260923426511168 -10.119661325070803 1.5927557243537434 -6.10662600432525 0.0 -2.4254235921784115 -1.9671309526934084 -3.3286180099265987 -2.8499659817055276 -1.7302731969932739 -1.1837682861777055 -0.11782983264763462 -0.3265440770990644 -1.0546003416785925 0.0 0.0 0.0 0.0; -1.6061109633450728 0.3963654501432292 -8.331365272058026 -5.264727157230805 -0.001 -7.627995401504769 -6.199926026836587 -10.51167941872737 -9.015447054009817 -5.5777954205960585 -3.778335116508615 -0.3202631508161873 -0.9869003787621273 -3.388660360774488 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 -3.548140511799198 -2.6482424961850004 -7.102082038003319 -6.208194955680467 0.0 0.0 0.0 0.0 0.0 -3.548140511799198 -2.6482424961850004 -8.102082038003319 -6.208194955680467; 0.0 0.0 0.0 0.0 0.0 -2.648242496185001 -2.0513640329992824 -5.597632716308486 -3.3299392371108265 0.0 0.0 0.0 0.0 0.0 -2.648242496185001 -2.0513640329992824 -5.597632716308486 -4.3299392371108265; 0.0 0.0 0.0 0.0 0.0 -18.102082038003317 4.4023672836915155 -24.873852981120933 -13.163260306077984 2.3828136105535203 2.0893223444460514 2.283219179019077 2.191129612074663 0.9999995532335706 -8.102082038003319 -5.5976327163084845 -21.873852981120933 -16.163260306077984; 0.0 0.0 0.0 0.0 0.0 3.7918050443195312 -14.329939237110827 -13.16326030607799 -15.408656206873678 0.0 0.0 0.0 0.0 0.0 -6.208194955680469 -4.3299392371108265 -16.16326030607799 -12.408656206873678; 0.0 0.0 0.0 0.0 0.0 -2.350321265186566 -1.6840654718182213 -5.594310418261385 -4.101569094377303 -1.488412883057869 -0.9639113397259123 0.5882419313046784 -0.8312794084435436 -0.3659595949709856 -1.5518179441443989 -1.040326764393575 -4.494409235292025 -3.159979998914917; 0.0 0.0 0.0 0.0 0.0 -2.2279953968524775 -1.5507658631630432 -5.630031465101207 -4.120845772910192 -0.9757754051184105 -1.0966951955892246 -0.4969559273832133 0.8240379338870026 -0.24851140114082376 -1.6974243978618049 -1.1252972716221998 -4.898931889735011 -3.4955631438299775; 0.0 0.0 0.0 0.0 0.0 -3.8534645754285197 -3.020606376036092 -6.73680504859119 -5.278229500140075 -14.905661141735358 5.981263475083093 -9.010086642052235 -0.969073690844096 -1.2936529189657013 -0.9518255041454229 -0.6650001476111296 -2.7386129367302487 -1.850027379296099; 0.0 0.0 0.0 0.0 0.0 -4.441954983917959 -3.3488668822504266 -9.320785111635601 -6.79570619414841 7.530650229517894 -11.303429611248509 1.4749258917061088 -6.433170081424314 -1.0546003416785925 -2.016531391739547 -1.381735929557018 -5.9921671017090015 -3.945740212442883; 0.0 0.0 0.0 0.0 0.0 -24.18765275168101 -17.451102917064617 -57.03910030189987 -42.6620899621136 -7.183906383941132 -3.3819696663653858 -8.651628422874213 -6.251627535992933 -3.3896603607744877 -16.55965735017624 -11.251176890228031 -46.5274208831725 -33.64664290810378; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -10.0 10.0 -3.0 3.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 10.0 -10.0 3.0 -3.0]
     marginplottestB = [0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 1.0; 0.0;;]
     marginplottestC = [2.3828136105535203 2.0893223444460514 2.283219179019077 2.191129612074663 0.9999995532335706 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]
     marginplottestD = [0.0;;]
     ss(marginplottestA, marginplottestB, marginplottestC, marginplottestD)
end

marginplot(G)
ωgₘ, gₘ, ωϕₘ, ϕₘ = margin(G, allMargins=true)
@test ωgₘ[] ≈ [0.9788420039794175, 21.96146912516854]
@test gₘ[] ≈ [0.19398511072183675, 11.886085245062574]
@test ωϕₘ[][] ≈ 3.484166418318219
@test ϕₘ[][] ≈ 50.783187269018754

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
@test linfnorm(gof[1] - sensitivity(P,C))[1] < 1e-9
@test linfnorm(gof[2] - G_PS(P,C))[1] < 1e-9
@test linfnorm(gof[3] - G_CS(P,C))[1] < 1e-9
@test linfnorm(gof[4] - comp_sensitivity(P,C))[1] < 1e-9
@test_nowarn gangoffourplot(P, C)
@test_nowarn gangoffourplot([P, P], C)
@test_nowarn gangoffourplot(P, [C, C])

F = ssrand(2,2,2,proper=true)
@test_nowarn gangofseven(P, C, F)

## Test bound functions
M_S = 1.3
g_m, ϕ_m = margin_bounds(M_S)
@test Ms_from_gain_margin(g_m) ≈ M_S
@test Ms_from_phase_margin(ϕ_m) ≈ M_S

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



## large TF poles and zeros
G = ssrand(2,3,4)
Gtf = tf(G)

pss  = poles(G)
zss  = tzeros(G)
ptf  = poles(Gtf)
ztf  = tzeros(Gtf)
pzpk = poles(zpk(G))
zzpk = tzeros(zpk(G))