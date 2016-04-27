module TestAnalysis
using CustomTest
using ControlSystems

## TZERO ##
# Examples from the Emami-Naeini & Van Dooren Paper
# Example 3
A = [0 1 0 0 0 0;
     0 0 1 0 0 0;
     0 0 0 0 0 0;
     0 0 0 0 1 0;
     0 0 0 0 0 1;
     0 0 0 0 0 0]
B = [0 0 1 0 0 0;
     0 0 0 0 0 1]'
C = [1 1 0 0 0 0;
     0 0 0 1 -1 0]
D = [1 0;
     1 0]

ex_3 = ss(A, B, C, D)
@test_approx_eq tzero(ex_3) [0.3411639019140099 + 1.161541399997252im,
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
@test_approx_eq tzero(ex_4) [-0.06467751189940692,-0.3680512036036696]

# Example 5
s = tf("s")
ex_5 = 1/s^15
@test tzero(ex_5) == Float64[]
@test tzero(ss(ex_5)) == Float64[]

# Example 6
A = [2 -1 0;
     0 0 0;
     -1 0 0]
B = [0 0 1]'
C = [0 -1 0]
D = [0]
ex_6 = ss(A, B, C, D)
@test tzero(ex_6) == Float64[]

# Example 7
ex_7 = ss(zeros(2, 2), [0 1]', [-1 0], [0])
@test tzero(ex_7) == Float64[]

# Example 8
A = [-2 1 0 0 0 0;
     1 -2 1 0 1 -1;
     0 1 -2 1 0 0;
     0 0 1 -1 0 1;
     0 -1 0 0 0 0;
     0 1 0 -1 0 0]
B = [1 0 0 0 1 0]'
C = [0 0 0 1 0 0]
D = [0]
ex_8 = ss(A, B, C, D)
# TODO : there may be a way to improve the precision of this example.
@test_approx_eq_eps tzero(ex_8) [-1.0, -1.0] 1e-7

# Example 9
ex_9 = (s - 20)/s^15
@test_approx_eq tzero(ex_9) [20.0]
@test_approx_eq tzero(ss(ex_9)) [20.0]

# Example 11
A = [-2 -6 3 -7 6;
     0 -5 4 -4 8;
     0 2 0 2 -2;
     0 6 -3 5 -6;
     0 -2 2 -2 5]
B = [-2 -8 -3 1 -8;
     7 -5 0 5 0]'
C = [0 -1 2 -1 -1;
     1 1 1 0 -1;
     0 3 -2 3 -1]
D = [0 0;
     0 0;
     0 0]
ex_11 = ss(A, B, C, D)
@test_approx_eq tzero(ex_11) [4.0, -3.0]

# Test for multiple zeros, siso tf
sys = s*(s + 1)*(s^2 + 1)*(s - 3)/((s + 1)*(s + 4)*(s - 4))
@test_approx_eq tzero(sys) [-1.0, -im, im, 3.0, 0.0]

## POLE ##
@test_approx_eq pole(sys) [-1.0, 4.0, -4.0]
@test_approx_eq pole([sys sys]) [-1.0, 4.0, -4.0, -1.0, 4.0, -4.0]
@test_approx_eq pole(ex_11) eig(ex_11.A)[1]

## ZPKDATA ##
# Sort a complex vector by real, breaking ties with imag
sortcomplex(a) = sort!(sort(a, by=imag), alg=MergeSort, by=real)
# Compare each vector in an array of vectors
macro test_array_vecs_eps(a, b, tol)
    quote
        @test size($a) == size($b)
        for (res, sol) = zip($a, $b)
            @test_approx_eq_eps sortcomplex(res) sol $tol
        end
    end
end
H = [tf(0) tf([3, 0],[1, 1, 10]) ; tf([1, 1],[1, 5]) tf([2],[1, 6])]
G = ss(H)
sol_z = vecarray(Complex128, 2, 2, Complex128[], Complex128[0.0 + 0.0im],
        Complex128[-1.0 + 0.0im], Complex128[])
sol_p = vecarray(Complex128, 2, 2, Complex128[], Complex128[-0.5 - 3.1224989991991996im,
        -0.5 + 3.1224989991991996im],
        Complex128[-5.0 + 0.0im], Complex128[-6.0 + 0.0im])
sol_k = [0.0 3.0; 1.0 2.0]
z, p, k = zpkdata(H)
@test_array_vecs_eps z sol_z 2*eps(Complex128)
@test_array_vecs_eps p sol_p 2*eps(Complex128)
@test k == sol_k
z, p, k = zpkdata(G)
@test_array_vecs_eps z sol_z 10*eps(Complex128)
@test_array_vecs_eps p sol_p 10*eps(Complex128)
@test k == sol_k

## GAIN ## #Gain is confusing when referring to zpkdata. Test dcgain instead
@test [dcgain(H[1, 1]) dcgain(H[1, 2]); dcgain(H[2, 1]) dcgain(H[2, 2])] ≈ [0 0; 0.2 1/3]
@test [dcgain(G[1, 1]) dcgain(G[1, 2]); dcgain(G[2, 1]) dcgain(G[2, 2])] ≈ [0 0; 0.2 1/3]
@test_err dcgain(G)

## MARKOVPARAM ##
@test markovparam(G, 0) == [0.0 0.0; 1.0 0.0]
@test markovparam(G, 1) == [0.0 3.0; -4.0 2.0]
@test markovparam(G, 2) == [0.0 -3.0; 20.0 -12.0]

## DAMP ##
@test_approx_eq damp(sys)[1] [1.0, 4.0, 4.0]
@test_approx_eq damp(sys)[2] [1.0, -1.0, 1.0]
@test_approx_eq damp(ex_11)[1] [1.0, 1.0, 2.0, 2.0, 3.0]
@test_approx_eq damp(ex_11)[2] [1.0, -1.0, -1.0, 1.0, -1.0]

## DAMPREPORT ##
@test sprint(dampreport, sys) == (
"|     Pole      |   Damping     |   Frequency   | Time Constant |\n"*
"|               |    Ratio      |   (rad/sec)   |     (sec)     |\n"*
"+---------------+---------------+---------------+---------------+\n"*
"|  -1.000e+00   |  1.000e+00    |  1.000e+00    |  1.000e+00    |\n"*
"|  4.000e+00    |  -1.000e+00   |  4.000e+00    |  -2.500e-01   |\n"*
"|  -4.000e+00   |  1.000e+00    |  4.000e+00    |  2.500e-01    |\n")
@test sprint(dampreport, ex_11) == (
"|     Pole      |   Damping     |   Frequency   | Time Constant |\n"*
"|               |    Ratio      |   (rad/sec)   |     (sec)     |\n"*
"+---------------+---------------+---------------+---------------+\n"*
"|  -1.000e+00   |  1.000e+00    |  1.000e+00    |  1.000e+00    |\n"*
"|  1.000e+00    |  -1.000e+00   |  1.000e+00    |  -1.000e+00   |\n"*
"|  2.000e+00    |  -1.000e+00   |  2.000e+00    |  -5.000e-01   |\n"*
"|  -2.000e+00   |  1.000e+00    |  2.000e+00    |  5.000e-01    |\n"*
"|  3.000e+00    |  -1.000e+00   |  3.000e+00    |  -3.333e-01   |\n")

end
