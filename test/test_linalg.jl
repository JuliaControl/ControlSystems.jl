@testset "test_linalg" begin
a = [-3 2;1 1]
b = [0  1]'
c = [1 -1]
r = 3

@test care(a, b, c'*c, r) ≈ [0.5895174372762604 1.8215747248860816; 1.8215747248860823 8.818839806923107]
@test dare(a, b, c'*c, r) ≈ [240.73344504496302 -131.09928700089387; -131.0992870008943 75.93413176750603]

## Test dare method for non invertible matrices
A = [0 1; 0 0];
B0 = [0; 1];
Q = Matrix{Float64}(I, 2, 2);
R0 = 1
X = dare(A, B0, Q, R0);
# Reshape for matrix expression
B = reshape(B0, 2, 1)
R = fill(R0, 1, 1)
@test norm(A'X*A - X - (A'X*B)*((B'X*B + R)\(B'X*A)) + Q) < 1e-14
## Test dare for scalars
A = 1.0;
B = 1.0;
Q = 1.0;
R = 1.0;
X0 = dare(A, B, Q, R);
X = X0[1]
@test norm(A'X*A - X - (A'X*B)*((B'X*B + R)\(B'X*A)) + Q) < 1e-14
# Test for complex matrices
I2 = Matrix{Float64}(I, 2, 2)
@test dare([1.0 im; im 1.0], I2, I2, I2) ≈ [1 + sqrt(2) 0; 0 1 + sqrt(2)]
# And complex scalars
@test dare(1.0, 1, 1, 1) ≈ fill((1 + sqrt(5))/2, 1, 1)
@test dare(1.0im, 1, 1, 1) ≈ fill((1 + sqrt(5))/2, 1, 1)
@test dare(1.0, 1im, 1, 1) ≈ fill((1 + sqrt(5))/2, 1, 1)

## Test gram, ctrb, obsv
a_2 = [-5 -3; 2 -9]
C_212 = ss(a_2, [1; 2], [1 0; 0 1], [0; 0])
C_222 = ss(a_2, [1 0; 0 2], [1 0; 0 1], zeros(2,2))

a_3 = [-0.81  0.47 -0.43  1.58  0.25 -0.40  0.92;
       -0.61 -1.89  0.80 -1.59  2.01  0.98 -0.90;
        0.50 -1.18 -2.12 -1.56 -1.14  0.14 -0.87;
       -1.28  2.12  0.47 -1.16  3.65 -1.23 -1.26;
       -0.24 -0.08  1.58 -3.58 -1.26  1.72 -2.59;
        1.28 -0.96 -1.28 -0.57 -2.43 -2.43 -0.36;
       -0.16  1.53 -0.99  1.47  0.61 -2.22 -3.27]
b_3 =  [ 0.49 -0.80  0.10;
        -0.18  0.70  0.72;
         0.00  0.83  2.58;
         1.42  0.00  0.00;
         0.00  0.22  0.19;
         0.00  0.00 -0.08;
         0.00 -1.15 -1.93]
c_3 = [ 0.00  0.84  0.10  0.30  0.49  1.71  0.00
       -1.79 -0.89 -0.54 -0.60  0.74 -0.19  0.00]
d_3 = [ 1.36  0.96  1.44
       -1.07  0.00 -1.96]
C_732 = ss(a_3,b_3,c_3,d_3);

s = tf("s")
f_C_211 = (s+2)*(s+3)/((s+4)*(s+5))
# biquad passband
omega0 = 52.0; Q = 10
f_C_211_bis = (s/(Q*omega0)) / ((s/omega0)^2 + s/(Q*omega0) + 1)    # Works again in julia 1.0

C_22tf = [0 tf([3,0],[1,1,10]);tf([1,1],[1,5]) tf(2,[1,6])]

da_2 = [0.2 -0.8; -0.8 0.07]
D_221 = ss(da_2, [1 0; 0 2], [1 0], [0 0], 0.005)
D_222 = ss(da_2, [1 0; 0 2], [1 0; 0 1], zeros(2,2), 0.005)

a_4 = [-0.08 0.83 0.0 0.0;-0.83 -0.08 0.0 0.0;0.0 0.0 -0.7 9.0;0.0 0.0 -9.0 -0.7]
b_4 = [1 1;0 0;1 -1;0 0]
c_4 = [0.4 0.0 0.4 0.0;0.6 0.0 1.0 0.0]
d_4 = [0.3 0.0;0.0 -0.15]
D_422 = ss(a_4,b_4,c_4,d_4,0.2)
z = tf("z", 0.05)
D_311 = z^3/((z+0.5)*(z^2-1.4z+0.72801))

C_static = ss([-0.78593 -1.2707;  0.2638 -1.13143])
D_static = ss([0.704473 1.56483; -1.6661 -2.16041], 0.07)


@test gram(C_212,:c) ≈ [0.042016806722689065 0.09663865546218485
    0.09663865546218488 0.24369747899159663]
@test gram(C_212,:o) ≈ [0.09523809523809523 -0.0119047619047619
    -0.011904761904761911 0.059523809523809534]
@test gram(C_222,:c) ≈ [0.11764705882352941 -0.029411764705882346
    -0.02941176470588234 0.21568627450980393]
@test gram(C_222,:o) ≈ [0.09523809523809523 -0.0119047619047619
    -0.011904761904761911 0.059523809523809534]
@test gram(D_221,:c) ≈ [11.895658388619264 -7.526846649259635
    -7.526846649259634 12.517564258299078]
@test gram(D_221,:o) ≈ [3.1039794541519283 -1.7910988952221152
    -1.7910988952221152 2.197919733616834]
@test gram(D_222,:c) ≈ [11.895658388619264 -7.526846649259635
    -7.526846649259634 12.517564258299078]
@test gram(D_222,:o) ≈ [5.301899187768763 -3.2250358337314955
    -3.2250358337314955 4.777830864787395]

@test obsv(C_222) == [1 0; 0 1; -5 -3; 2 -9]
@test ctrb(C_222) == [1 0 -5 -6; 0 2 2 -18]
@test obsv(C_212) == [1 0; 0 1; -5 -3; 2 -9]
@test ctrb(C_212) == [1 -11; 2 -16]
@test norm(C_222) ≈ 0.5773502691896258
@test norm(C_212) ≈ 0.5345224838248488
@test norm(D_222) ≈ 4.940973856125768
@test norm(D_221) ≈ 3.4490083195926426
@test norm(ss([1],[2],[3],[4])) == Inf

#
## Test Hinfinity norm computations
#
# Continuous
Ninf, ω_peak = hinfnorm(ss(-1, 1, 1, 0))
@test Ninf ≈ 1 rtol=1e-6
@test ω_peak ≈ 0 rtol=1e-6

Ninf, ω_peak = hinfnorm(ss(-1, 1, 1, 1))
@test Ninf ≈ 2 rtol=1e-6
@test ω_peak ≈ 0 rtol=1e-6

Ninf, ω_peak = hinfnorm(ss(1.5))
@test Ninf ≈ 1.5 rtol=1e-6
@test ω_peak ≈ 0 rtol=1e-6


tolHinf = 1e-12
@test norm(C_212, Inf, tol=tolHinf) ≈ 0.242535625036333 atol=tolHinf
@test norm(C_222, Inf, tol=tolHinf) ≈ 0.242535625036333 atol=tolHinf
Ninf, ω_peak = hinfnorm(C_732, tol=tolHinf)
@test Ninf ≈ 4.899135403568278 atol=(10*tolHinf)
@test ω_peak ≈ 6.112977387441163 atol=1e-6
@test norm(f_C_211, Inf, tol=tolHinf) ≈ 1.0 atol=(2*tolHinf)
Ninf, ω_peak =  hinfnorm(f_C_211_bis, tol=tolHinf)
@test Ninf ≈ 1.0
@test ω_peak ≈ 52.0

# unstable system
sys = 1/(s-1)
# @inferred hinfnorm(sys)
@test norm(sys, Inf, tol=tolHinf) ≈ Inf
Ninf, ω_peak = hinfnorm(sys; tol=tolHinf)
@test Ninf == Inf
@test isnan(ω_peak)
Ninf, ω_peak = linfnorm(sys; tol=tolHinf)
@test Ninf ≈ 1
@test ω_peak ≈ 0


Ninf, ω_peak = hinfnorm(C_static, tol=tolHinf)
@test Ninf ≈ 1.760164138376307 atol=1e-12
@test ω_peak ≈ 0 atol=1e-12
# @inferred hinfnorm(C_static, tol=tolHinf)

Ninf, ω_peak = hinfnorm(C_22tf, tol=tolHinf)
@test Ninf ≈ 3.014974550173459 atol=(10*tolHinf)
@test ω_peak ≈ 3.162123338668049 atol=1e-8

# Poles on the imaginary axis
# Only integrator
Ninf, ω_peak = hinfnorm(ss(0.0, 1, 1, 0))
@test Ninf == Inf
@test ω_peak ≈ 0

# Integrator and lag
Ninf, ω_peak = hinfnorm(ss([-2.0 1; 0 0], [0; 1], [1 0], 0))
@test Ninf == Inf
@test ω_peak ≈ 0

# 1/(s^2 + 1)
Ninf, ω_peak = hinfnorm(ss([0 1.0; -1.0 0], [1; 0], [0 1], 0))
@test Ninf == Inf
@test ω_peak ≈ 1

# 1/(s - im), complex coefficients
Ninf, ω_peak = hinfnorm(ss(im, 1, 1, 0))
@test Ninf == Inf
@test ω_peak ≈ 1


# Complex coefficient systems
Ninf, ω_peak = hinfnorm(ss(-1+im, 1, 1, 0))
@test Ninf ≈ 1 rtol=1e-6
@test ω_peak ≈ 1 rtol=1e-3

Ninf, ω_peak = hinfnorm(ss(-1-im, 1, 1, 0))
@test Ninf ≈ 1 rtol=1e-6
@test ω_peak ≈ -1 rtol=1e-3

# 1/(s + 1)/(s - 1)
@test hinfnorm(ss([-1 1; 0 1], [0; 1], [1 0], 0))[1] == Inf
@test linfnorm(ss([-1 1; 0 1], [0; 1], [1 0], 0))[1] ≈ 1 rtol=1e-6

# Second order resonant system with peaks at ±sqrt(2)
G = ss([0 1; -1 -1], [0; 1], [1 0], 0)
Ninf, ω_peak = hinfnorm(G)
@test Ninf ≈ 2/sqrt(3) rtol=1e-6
@test ω_peak ≈ 1/sqrt(2) rtol=1e-3
@test hinfnorm(G, tol=1e-8)[1] ≈ 2/sqrt(3) rtol=1e-8
@test hinfnorm(G, tol=1e-10)[1] ≈ 2/sqrt(3) rtol=1e-10

# Same version as above but frequency shifted by -1, i.e., two peaks at -1 ± sqrt(2)
Ninf, ω_peak = hinfnorm(ss([0 1; -im -(1+2im)], [0; 1], [1 0], 0))
@test Ninf ≈ 2/sqrt(3) rtol=1e-6
@test abs(ω_peak+1) ≈ 1/sqrt(2) rtol=1e-3

# Same version as above but frequency shifted by +1, i.e., two peaks at 1 ± sqrt(2)
Ninf, ω_peak = hinfnorm( ss([0 1; im -1+2im], [0; 1], [1 0], 0))
@test Ninf ≈ 2/sqrt(3) rtol=1e-6
@test abs(ω_peak-1) ≈ 1/sqrt(2) rtol=1e-3


# System with mutiple resonance peaks
A = [-0.1*I        -Diagonal(1:10)
     Diagonal(1:10) zeros(10,10)]

Ninf, ω_peak = hinfnorm(ss(A, 0.1*[ones(10); zeros(10)], [zeros(1,10) ones(1,10)], 0))
@test Ninf ≈ 1.05871265 rtol=1e-5
@test ω_peak ≈ 0.9874 rtol=1e-4

Ninf, ω_peak = hinfnorm(ss(A, 0.1*[ones(4); 2; ones(5); zeros(10)], [zeros(1,10) ones(1,10)], 0))
@test Ninf ≈ 2.00918 rtol=1e-5
@test ω_peak ≈ 4.9981 rtol=1e-4



# Discrete Time
Ninf, ω_peak =  hinfnorm(ss(0.8, 1, 1, 0, 1))
@test Ninf ≈ 5 rtol=1e-6
@test ω_peak ≈ 0

Ninf, ω_peak =  hinfnorm(ss(-0.8, 1, 1, 0, 1))
@test Ninf ≈ 5 rtol=1e-6
@test ω_peak ≈ pi

Ninf, ω_peak = hinfnorm(D_221, tol=tolHinf)
@test Ninf ≈ 17.794697451669421 atol=(20*tolHinf)
@test ω_peak ≈ 0 atol=1e-8

Ninf, ω_peak = linfnorm(D_422, tol=tolHinf)
@test Ninf ≈ 3.360351099392252 atol=(10*tolHinf)
@test ω_peak ≈ 8.320643111730551 atol=1e-6
Ninf, ω_peak = hinfnorm(D_422, tol=tolHinf)
@test Ninf == Inf
@test isnan(ω_peak)


Ninf, ω_peak = hinfnorm(D_311, tol=tolHinf)
@test Ninf ≈ 4.458729529942810 atol=(10*tolHinf)
@test ω_peak ≈ 11.878021287349698 atol=1e-6

Ninf, ω_peak = hinfnorm(D_static, tol=tolHinf)
@test Ninf ≈ 3.205246234285972 atol=1e-12
@test ω_peak ≈ 0 atol=1e-12

# Pole on unit circle
Ninf, ω_peak = hinfnorm(ss(1, 1, 1, 0, 1))
@test Ninf == Inf
@test ω_peak ≈ 0

# Complex coefficients
Ninf, ω_peak =  hinfnorm(ss(0.8*exp(1im), 1, 1, 0, 1))
@test Ninf ≈ 5 rtol=1e-6
@test ω_peak ≈ 1

Ninf, ω_peak =  hinfnorm(ss(0.8*exp(-1im), 1, 1, 0, 1))
@test Ninf ≈ 5 rtol=1e-6
@test ω_peak ≈ -1

sys = hinfnorm(ss(0.8*exp(2.5im), 1, 1, 0, 1))
Ninf, ω_peak =  hinfnorm(ss(0.8*exp(2.5im), 1, 1, 0, 1))
@test Ninf ≈ 5 rtol=1e-6
@test ω_peak ≈ 2.5
Ninf, ω_peak = linfnorm(ss(0.8*exp(2.5im), 1, 1, 0, 1))
@test Ninf ≈ 5 rtol=1e-6
@test ω_peak ≈ 2.5


# Unstable discrete-time system
sys = ss(2, 3, 1, 1, 1)
Ninf, ω_peak = hinfnorm(sys)
@test Ninf == Inf
@test isnan(ω_peak)
Ninf, ω_peak = linfnorm(sys)
@test Ninf ≈ 2.0 rtol=1e-6
@test ω_peak ≈ 0


A = [1  100  10000; .01  1  100; .0001  .01  1]
T, P, B = balance(A)
# The scaling is BLAS dependent. However, the ratio should be the same on all
# machines. We just need to check that T == res * constant
res_diag = [512, 8, 0.0625]
constant = T[1]/res_diag[1]
@test all(diag(T) == res_diag * constant)
@test P == [1 0 0; 0 1 0; 0 0 1]
@test B ≈ [1.0 1.5625 1.220703125; 0.64 1.0 0.78125; 0.8192 1.28 1.0]
end
