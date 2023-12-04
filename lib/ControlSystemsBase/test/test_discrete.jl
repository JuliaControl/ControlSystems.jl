@testset "test_discrete" begin

C_111 = ss([-5], [2], [3], [0])
C_212 = ss([-5 -3; 2 -9], [1; 2], [1 0; 0 1], [0; 0])
C_221 = ss([-5 -3; 2 -9], [1 0; 0 2], [1 0], [0 0])
C_222_d = ss([-5 -3; 2 -9], [1 0; 0 2], [1 0; 0 1], [1 0; 0 1])

@test c2d_x0map(ss(4*[1 0; 0 1]), 0.5, :zoh) == (ss(4*[1 0; 0 1], 0.5), zeros(0, 2))
@test c2d_x0map(ss(4*[1 0; 0 1]), 0.5, :foh) == (ss(4*[1 0; 0 1], 0.5), zeros(0, 2))
@test_tupleapprox(c2d_x0map(C_111, 0.01, :zoh),
    ss([0.951229424500714], [0.019508230199714396], [3], [0], 0.01), [1 0])
@test_tupleapprox(c2d_x0map(C_111, 0.01, :foh),
    ss([0.951229424500714], [0.01902855227625244], [3], [0.029506188017136226], 0.01),
    [1 -0.009835396005712075])
@test_tupleapprox(c2d_x0map(C_212, 0.01, :zoh),
    ss([0.9509478368863918 -0.027970882212682433; 0.018647254808454958 0.9136533272694819],
    [0.009466805409666932; 0.019219966830212765], [1 0; 0 1], [0; 0], 0.01),
    [1 0 0; 0 1 0])
@test_tupleapprox(c2d_x0map(C_212, 0.01, :foh),
    ss([0.9509478368863921 -0.027970882212682433; 0.018647254808454954 0.913653327269482],
    [0.008957940478201584; 0.018468989584974498], [1 0; 0 1], [0.004820885889482196;
    0.009738343195298675], 0.01), [1 0 -0.004820885889482196; 0 1 -0.009738343195298675])
@test_tupleapprox(c2d_x0map(C_221, 0.01, :zoh),
    ss([0.9509478368863918 -0.027970882212682433; 0.018647254808454958 0.9136533272694819],
    [0.009753161420545834 -0.0002863560108789034; 9.54520036263011e-5 0.019124514826586465],
    [1.0 0.0], [0.0 0.0], 0.01), [1 0 0 0; 0 1 0 0])
@test_tupleapprox(c2d_x0map(C_221, 0.01, :foh),
    ss([0.9509478368863921 -0.027970882212682433; 0.018647254808454954 0.913653327269482],
    [0.009511049106772921 -0.0005531086285713394; 0.00018436954285711309 0.018284620042117387],
    [1 0], [0.004917457305816479 -9.657141633428213e-5], 0.01),
    [1 0 -0.004917457305816479 9.657141633428213e-5;
    0 1 -3.219047211142736e-5 -0.009706152723187249])
@test_tupleapprox(c2d_x0map(C_222_d, 0.01, :zoh),
    ss([0.9509478368863918 -0.027970882212682433; 0.018647254808454958 0.9136533272694819],
    [0.009753161420545834 -0.0002863560108789034; 9.54520036263011e-5 0.019124514826586465],
    [1 0; 0 1], [1 0; 0 1], 0.01), [1 0 0 0; 0 1 0 0])
@test_tupleapprox(c2d_x0map(C_222_d, 0.01, :foh),
    ss([0.9509478368863921 -0.027970882212682433; 0.018647254808454954 0.913653327269482],
    [0.009511049106772921 -0.0005531086285713394; 0.00018436954285711309 0.018284620042117387],
    [1 0; 0 1], [1.0049174573058164 -9.657141633428213e-5; 3.219047211142736e-5 1.0097061527231872], 0.01),
    [1 0 -0.004917457305816479 9.657141633428213e-5; 0 1 -3.219047211142736e-5 -0.009706152723187249])

# Test c2d for transfer functions
G = tf([1, 1], [1, 3, 1])
Gd = c2d(G, 0.2)
@test typeof(Gd) <: TransferFunction{<:Discrete, <:ControlSystemsBase.SisoRational}
@test Gd ≈ tf([0, 0.165883310712090, -0.135903621603238], [1.0, -1.518831946985175, 0.548811636094027], 0.2) rtol=1e-14

# Issue #391
@test c2d(tf(C_111), 0.01, :zoh) ≈ tf(c2d(C_111, 0.01, :zoh))
@test c2d(tf(C_111), 0.01, :foh) ≈ tf(c2d(C_111, 0.01, :foh))

# Issue #471 discussion, test c2d for zero row C
C_210 = ss(C_212.A, C_212.B, zeros(0, 2), zeros(0, 1))
@test c2d(C_210, 0.01).A ≈ c2d(C_212, 0.01).A

# c2d on a zpk model should arguably return a zpk model
Gdzpk = c2d(zpk(G), 0.2)
@test typeof(Gdzpk) <: TransferFunction{<:Discrete, <:ControlSystemsBase.SisoZpk}
@test Gdzpk ≈ zpk([0.8192724211968178], [0.9264518519788194, 0.5923800950063556], 0.16588331071209, 0.2) rtol=1e-14


# ERRORS
@test_throws MethodError c2d(ss([1], [2], [3], [4], 0.01), 0.01)   # Already discrete
@test_throws MethodError c2d(ss([1], [2], [3], [4], -1), 0.01)     # Already discrete


# d2c
@test d2c(c2d(C_111, 0.1)) ≈ C_111
@test d2c(c2d(C_212, 0.1)) ≈ C_212
@test d2c(c2d(C_221, 0.1)) ≈ C_221
@test d2c(c2d(C_222_d, 0.1)) ≈ C_222_d
@test d2c(Gd) ≈ G

sys = ss([0 1; 0 0], [0;1], [1 0], 0)
sysd = c2d(sys, 1)
@test d2c(sysd) ≈ sys


# forward euler / tustin
@test c2d(C_111, 1, :fwdeuler).A == I + C_111.A
for method in (:fwdeuler, :tustin)
    @test d2c(c2d(C_111, 0.01, method), method) ≈ C_111 atol = sqrt(eps())
    @test d2c(c2d(C_212, 0.01, method), method) ≈ C_212 atol = sqrt(eps())
    @test d2c(c2d(C_221, 0.01, method), method) ≈ C_221 atol = sqrt(eps())
    @test d2c(c2d(C_222_d, 0.01, method), method) ≈ C_222_d atol = sqrt(eps())
    @test d2c(c2d(G, 0.01, method), method) ≈ G atol = sqrt(eps())
end

matlab_tustin = let 
    A = [0.95094630230333 -0.02800401390866; 0.018669342605773 0.913607617091783]
    B = [0.009754731511517 -0.000280040139087; 9.3346713029e-5 0.019136076170918]
    C = [0.975473151151665 -0.01400200695433; 0.009334671302887 0.956803808545892]
    D = [1.004877365755758 -0.000140020069543; 4.6673356514e-5 1.009568038085459]
    ss(A,B,C,D,0.01)
end

sys_tustin = c2d(C_222_d, 0.01, :tustin)
@test sys_tustin ≈ matlab_tustin atol = 1e-12

matlab_prewarp = let
    A = [0.950906169509825 -0.028025790680969; 0.018683860453979 0.913538448601866]
    B = [0.009762667760265 -0.00028049168885; 9.3497229617e-5 0.019151346602063]
    C = [0.975453084754912 -0.014012895340485; 0.00934193022699 0.956769224300933]
    D = [1.004881333880133 -0.000140245844425; 4.6748614808e-5 1.009575673301031]
    sys = ss(A,B,C,D,0.01)
end

sys_prewarp = c2d(C_222_d, 0.01, :tustin, w_prewarp=10)
@test sys_prewarp ≈ matlab_prewarp atol = 1e-12




# dab
import DSP: conv
ζ = 0.2
ω = 1
B = [1]
A   = [1, 2ζ*ω, ω^2]
P  = tf(B,A)
# Control design
ζ0 = 0.7
ω0 = 2
Am = [1, 2ζ0*ω0, ω0^2]
Ao = conv(2Am, [1/2, 1]) # Observer polynomial, add extra pole due to the integrator
AR = [1,0] # Force the controller to contain an integrator ( 1/(s+0) )

B⁺  = [1] # The process numerator polynomial can be facored as B = B⁺B⁻ where B⁻ contains the zeros we do not want to cancel (non-minimum phase and poorly damped zeros)
B⁻  = [1]
Bm  = conv(B⁺, B⁻) # In this case, keep the entire numerator polynomial of the process

R,S,T = rstc(B⁺,B⁻,A,Bm,Am,Ao,AR) # Calculate the 2-DOF controller polynomials

Gcl = tf(conv(B,T),zpconv(A,R,B,S)) # Form the closed loop polynomial from reference to output, the closed-loop characteristic polynomial is AR + BS, the function zpconv takes care of the polynomial multiplication and makes sure the coefficient vectors are of equal length

@test ControlSystemsBase.isstable(Gcl)

p = poles(Gcl)
# Test that all desired poles are in the closed-loop system
@test norm(minimum(abs.((poles(tf(Bm,Am)) .- sort(p, by=imag)')), dims=2)) < 1e-6
# Test that the observer poles are in the closed-loop system
@test norm(minimum(abs.((poles(tf(1,Ao)) .- sort(p, by=imag)')), dims=2)) < 1e-6


pd = c2d_poly2poly(A, 0.1)
@test pd ≈ denvec(c2d(P, 0.1))[]


@testset "d2c_exact" begin
    @info "Testing d2c_exact"
    sys = ssrand(2, 3, 2, Ts = 1)
    sysc = d2c_exact(sys)

    # bodeplot(sys, w, lab="Discrete SS")
    # bodeplot!(sysc, w, lab="Continuous SS with delays", l=:dash, size=(800, 800))
    w = exp10.(LinRange(-2, 2, 200))
    @test freqresp(sys, w) ≈ freqresp(sysc, w) atol = 1e-8
end


## Sampling of covariance and cost matrices

# Covariance matrices (opt = :o)
sysc = DemoSystems.resonant()
Ts = 0.5
sysd = c2d(sysc, Ts)
Qd = [1 0.1; 0.1 2]
Qc = d2c(sysd, Qd)

Qd2 = c2d(sysc, Qc, Ts) # Continuous system as input
@test Qd2 ≈ Qd

Qd3 = c2d(sysd, Qc) # Discrete system as input
@test Qd3 ≈ Qd


# Test sampling of cost matrix (opt = :c)
sysc = DemoSystems.resonant()
x0 = ones(sysc.nx)
Ts = 0.01 # cost approximation becomes more crude as Ts increases
Qc = [1 0.01; 0.01 2]
Rc = I(1)
sysd = c2d(sysc, Ts)

Qd, Rd = c2d(sysc, Qc, Rc, Ts, opt=:c) # Continuous system as input
Qc2, Rc2 = d2c(sysd, Qd, Rd; opt=:c) # Test round trip
@test Qc2 ≈ Qc
@test Rc2 ≈ Rc

Qd3, Rd3 = c2d(sysd, Qc, Rc; opt=:c) # Discrete system as input
@test Qd3 ≈ Qd
@test Rd3 ≈ Rd


# NOTE: these tests are not run due to OrdinaryDiffEq latency, they should pass
# using OrdinaryDiffEq
# L = lqr(sysc, Qc, Rc)
# dynamics = function (xc, p, t)
#     x = xc[1:sysc.nx]
#     u = -L*x
#     dx = sysc.A*x + sysc.B*u
#     dc = dot(x, Qc, x) + dot(u, Rc, u)
#     return [dx; dc]
# end
# prob = ODEProblem(dynamics, [x0; 0], (0.0, 10.0))
# sol = solve(prob, Tsit5(), reltol=1e-8, abstol=1e-8)
# cc = sol.u[end][end]
# Ld = lqr(sysd, Qd, Rd)
# sold = lsim(sysd, (x, t) -> -Ld*x, 0:Ts:10, x0 = x0)
# function cost(x, u, Q, R)
#     dot(x, Q, x) + dot(u, R, u)
# end
# cd = cost(sold.x, sold.u, Qd, Rd)
# @test cc ≈ cd rtol=0.01
# @test abs(cc-cd) < 1.0001*0.005531389319983315


# test case from paper
A = [
    2 -8 -6
    10 -19 -12
    -10 15 8
]
B = [
    5 1
    1 4
    3 2
]
Qc = [
    4 1 2
    1 3 1
    2 1 5
]

Qd = c2d(ss(A,B,I,0), Qc, 1, opt=:c)
Qd_van_load = [
    9.934877720 -11.08568953 -9.123023900
    -11.08568953 13.66870748 11.50451512
    -9.123023900 11.50451512 10.29179555 
]

@test norm(Qd - Qd_van_load) < 1e-6


end