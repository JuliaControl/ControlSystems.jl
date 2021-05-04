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
@test Gd ≈ tf([0, 0.165883310712090, -0.135903621603238], [1.0, -1.518831946985175, 0.548811636094027], 0.2) rtol=1e-14

# Issue #391
@test c2d(tf(C_111), 0.01, :zoh) ≈ tf(c2d(C_111, 0.01, :zoh))
@test c2d(tf(C_111), 0.01, :foh) ≈ tf(c2d(C_111, 0.01, :foh))

# Issue #471 discussion, test c2d for zero row C
C_210 = ss(C_212.A, C_212.B, zeros(0, 2), zeros(0, 1))
@test c2d(C_210, 0.01).A ≈ c2d(C_212, 0.01).A

# c2d on a zpk model should arguably return a zpk model
@test_broken typeof(c2d(zpk(G), 1)) <: TransferFunction{<:ControlSystems.SisoZpk}



# ERRORS
@test_throws MethodError c2d(ss([1], [2], [3], [4], 0.01), 0.01)   # Already discrete
@test_throws MethodError c2d(ss([1], [2], [3], [4], -1), 0.01)     # Already discrete


# d2c
@static if VERSION > v"1.4" # log(matrix) is buggy on previous versions, should be fixed in 1.4 and back-ported to 1.0.6
    @test d2c(c2d(C_111, 0.01)) ≈ C_111
    @test d2c(c2d(C_212, 0.01)) ≈ C_212
    @test d2c(c2d(C_221, 0.01)) ≈ C_221
    @test d2c(c2d(C_222_d, 0.01)) ≈ C_222_d
    @test d2c(Gd) ≈ G

    sys = ss([0 1; 0 0], [0;1], [1 0], 0)
    sysd = c2d(sys, 1)
    @test d2c(sysd) ≈ sys
end

# forward euler
@test c2d(C_111, 1, :fwdeuler).A == I + C_111.A
@test d2c(c2d(C_111, 0.01, :fwdeuler), :fwdeuler) ≈ C_111
@test d2c(c2d(C_212, 0.01, :fwdeuler), :fwdeuler) ≈ C_212
@test d2c(c2d(C_221, 0.01, :fwdeuler), :fwdeuler) ≈ C_221
@test d2c(c2d(C_222_d, 0.01, :fwdeuler), :fwdeuler) ≈ C_222_d
@test d2c(c2d(G, 0.01, :fwdeuler), :fwdeuler) ≈ G


Cd = c2d(C_111, 0.001, :fwdeuler)
t = 0:0.001:2
y,_ = step(C_111, t)
yd,_ = step(Cd, t)
@test norm(yd-y) / norm(y) < 1e-3

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

Gcl = tf(conv(B,T),zpconv(A,R,B,S)) # Form the closed loop polynomial from reference to output, the closed-loop characteristic polynomial is AR + BS, the function zpconv takes care of the polynomial multiplication and makes sure the coefficient vectores are of equal length

@test ControlSystems.isstable(Gcl)

p = pole(Gcl)
# Test that all desired poles are in the closed-loop system
@test norm(minimum(abs.((pole(tf(Bm,Am)) .- sort(p, by=imag)')), dims=2)) < 1e-6
# Test that the observer poles are in the closed-loop system
@test norm(minimum(abs.((pole(tf(1,Ao)) .- sort(p, by=imag)')), dims=2)) < 1e-6


end
