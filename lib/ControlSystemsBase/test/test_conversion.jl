@testset "test_conversion" begin

A = randn(3,3)
@inferred ControlSystemsBase.charpoly(A)

A = randn(ComplexF64,3,3)
@inferred ControlSystemsBase.charpoly(A)

G = tf(1.0,[1,1])
H = zpk([0.0], [1.0], 1.0)
# @inferred ControlSystemsBase.siso_tf_to_ss(Float64, G.matrix[1,1])
# @inferred ControlSystemsBase.siso_tf_to_ss(Float64, H.matrix[1,1])

# Easy second order system
sys1 = ss([-1 0;1 1],[1;0],[1 1],0)
G1 = tf([1,0],[1,0,-1])
H1 = zpk([0.0], [-1, 1.0], 1.0)

w = 10.0 .^ (-2:2:50)

@test ss(sys1) == sys1
@test freqresp(ss(G1), w) ≈ freqresp(sys1, w) rtol=1e-15
@test freqresp(ss(H1), w) ≈ freqresp(sys1, w) rtol=1e-15

# Test range in freqresp
@test freqresp(ss(G1), collect(1:2:10)) == freqresp(ss(G1), collect(1:2:10))

@test tf(sys1) ≈ G1 rtol=1e-15
@test tf(G1) == G1
@test tf(H1) == G1

@test zpk(sys1) ≈ H1 rtol=1e-15
@test zpk(G1) == H1
@test zpk(H1) == H1

# With direct term second order system
sys2 = ss([-2.5 0;1 1.5],[1;3],[1 2],2.5)
G2 = tf([2.5, 9.5, 6.125],[1, 1, -3.75])
H2 = zpk([-2.97703296142690, -0.82296703857310], [1.5, -2.5], 2.5)

w = 10.0 .^ (-2:2:50)

@test ss(sys2) == sys2
@test freqresp(ss(G2), w) ≈ freqresp(sys2, w) rtol=1e-15
@test freqresp(ss(H2), w) ≈ freqresp(sys2, w) rtol=1e-15

@test tf(sys2) ≈ G2 rtol=1e-15
@test tf(G2) == G2
@test tf(H2) ≈ G2 rtol=1e-15

@test zpk(sys2) ≈ H2 rtol=1e-15
@test zpk(G2) ≈ H2 rtol=1e-15
@test zpk(H2) == H2

# Non-proper system
G3 = tf([4.0, 1, 5],[2, 3])
H3 = zpk([(-1+im*sqrt(79))/8, (-1-im*sqrt(79))/8], [-3/2], 2)

@test tf(G3) == G3
@test tf(H3) ≈ G3 rtol=1e-2

@test zpk(G3) ≈ H3 rtol=1e-15
@test zpk(H3) == H3

s = tf("s")
@test_throws ErrorException delay(5)*((s+1))/((s+2)*(s+0.5))
@test_throws ErrorException delay(5)*zpk((s+1))/zpk((s+2)*(s+0.5))

P = ss(-1, 1, 1, 0)
Pd = P * delay(big(1))
@test Pd.P.P ≈ (tf(P) * delay(big(1))).P.P
@test Pd.Tau == [1]

# Test complex 1
A = [1.0 + im 1; 0 -2-3im]
B = [0;2]
C = [1 0]
D = zeros(1,1);
sys = ss(A,B,C,D)
sys2 = convert(TransferFunction, sys)
w = 10.0 .^ (-2:2:50)
@test freqresp(sys, w) ≈ freqresp(sys2, w)
csort = v -> sort(v, lt = (x,y) -> abs(x) < abs(y))
@test csort(poles(zpk(sys2))) ≈ [1+im, -2.0-3im]

# Test complex 2
sys4 = ss([-1.0+im 0;1 1],[1;0],[1 im],im)
tf(sys4)
G4 = tf([1.0im,2,-2],[1,-1im,-1+im])
H4 = zpk([im*(1 + sqrt(1+2im)), im*(1 - sqrt(1+2im))], [-1+im,1], 1.0im)

w = 10.0 .^ (-2:2:50)

@test ss(sys4) == sys4
@test freqresp(ss(G4), w) ≈ freqresp(sys4, w) rtol=1e-15
@test freqresp(ss(H4), w) ≈ freqresp(sys4, w) rtol=1e-15

@test tf(sys4) ≈ G4 rtol=1e-15
@test tf(G4) == G4
@test tf(H4) ≈ G4 rtol=1e-15

@test zpk(sys4) ≈ H4 rtol=1e-15
@test zpk(G4) ≈ H4 rtol=1e-15
@test zpk(H4) == H4

# Test conversion of Int system with direct term
w = [0.5, 1, 1.5, 2.0]
G = tf([2, 3],[5, 1])
@test freqresp(ss(G), w) ≈ freqresp(G, w)


sys2 = ss([-1 0;1 1],[1;0],[1 1],1)
G2 = tf([1,1,-1],[1,0,-1])
sys3 = ss([0 1;1 0],[0;1],[1 1],0)
G3 = tf([1,1],[1,0,-1])

#Test conversion of StateSpace to TransferFunction
@test tf(sys2) ≈ G2 atol=1e-15
@test tf(sys3) ≈ G3 atol=1e-15

#Test TransferFunction -> StateSpace
#ss is not unique so go back and forth
@test G2 ≈ tf(ss(G2)) atol=1e-15
@test G3 ≈ tf(ss(G3)) atol=1e-15

### test mimo
@test hinfnorm(zpk([sys1 sys2])-zpk([ss(sys1) ss(sys2)]))[1] < 1e-12


## Test some BigFloat
SSBigFloat = StateSpace{Continuous,BigFloat}
ZpkBigFloat = TransferFunction{Continuous,ControlSystemsBase.SisoZpk{BigFloat,Complex{BigFloat}}}
RationalBigFloat = TransferFunction{Continuous,ControlSystemsBase.SisoRational{BigFloat}}

s = tf("s")
f = zpk(1.0*(2s+3)/((5s+7)*(11s+13)))
fb = BigFloat(1.0)*f

# Reasonable defaults
# @test fb isa RationalBigFloat # Do we want this? We get TransferFunction{ControlSystemsBase.SisoZpk{BigFloat,Complex{Float64}}}
@test ss(fb) isa SSBigFloat
@test tf(fb) isa RationalBigFloat
# Can't compute poles:
# @test zpk(tf(fb)) isa SSBigFloat

@test convert(SSBigFloat, fb) isa SSBigFloat
@test convert(ZpkBigFloat, fb) isa ZpkBigFloat
@test convert(RationalBigFloat, fb) isa RationalBigFloat

## Test some rather trivial conversions of numeric types

b = 1.5
D22 = [1.0 2.0; 3.0 4.0]

@test convert(StateSpace{Continuous,Float64}, D22) == ss(D22)
@test convert(StateSpace{Continuous,Float64}, b) == ss(b)

@test convert(TransferFunction{Continuous,ControlSystemsBase.SisoRational{Float64}}, D22) == tf(D22)
@test convert(TransferFunction{Continuous,ControlSystemsBase.SisoRational{Float64}}, b) == tf(b)

@test convert(TransferFunction{Continuous,ControlSystemsBase.SisoZpk{Float64,ComplexF64}}, D22) == zpk(D22)
@test convert(TransferFunction{Continuous,ControlSystemsBase.SisoZpk{Float64,ComplexF64}}, b) == zpk(b)

## Test some discrete conversions
h = 0.1
sysd = ss([0.5 1; 0 1], [0; 1], [1 0], 0, h)
Gd = tf([1], [1, -1.5, 0.5], h)
Hd = zpk([], [1, 0.5], 1.0, h)

@test tf(sysd) ≈ Gd rtol=1e-15
@test tf(Gd) == Gd
@test tf(Hd) == Gd

@test zpk(sysd) ≈ Hd rtol=1e-15
@test zpk(Gd) == Hd
@test zpk(Hd) == Hd


# Error, not proper
@test_throws ControlSystemsBase.ImproperException ss(tf([1,0],[1]))

s1 = HeteroStateSpace(zeros(Int, 1,1),zeros(Int, 1,1),zeros(Int, 1,1),zeros(Int, 1,1))
s2 = HeteroStateSpace(zeros(Float64, 1,1),zeros(Float64, 1,1),zeros(Float64, 1,1),zeros(Float64, 1,1))
s3,s4 = promote(s1,s2)
@test s3 == s4 == s2
@test typeof(s3) == typeof(s4) == typeof(s2)

# https://github.com/JuliaControl/ControlSystems.jl/issues/828
mytf = [tf(1,[1,1]) tf(0); tf(0) tf(1,[1,1])]
myss = ss(mytf)
@test myss == ss([-1 0; 0 -1],[1 0; 0 1],[1 0; 0 1],[0 0; 0 0])

test = ss([-1 2 3; 0 -2 3; 0 0 -3],[1;2;3],[1 2 3], [0])
@test tf(test) isa TransferFunction


# Issue https://github.com/JuliaControl/ControlSystems.jl/issues/869

sys = tf(ssrand(2,2,2))
sys.matrix[1,1] = 0
syszpk = zpk(sys)
@test syszpk isa TransferFunction{Continuous, ControlSystemsBase.SisoZpk{Float64, ComplexF64}}

## Test multiplication of improper transfer function and statespace system
P = DemoSystems.double_mass_model()
C = tf('s')
@test norm(P*C - C*P) < 1e-10
mp = minreal(minreal(tf(P)*C) - P*C)
@test hinfnorm(mp)[1] < 1e-10
@inferred P*C

@test_logs (:warn,"Possible numerical instability detected: Multiplication of a statespace system and a non-proper transfer function may result in numerical inaccuracy. Verify result carefully, and consider making use of DescriptorSystems.jl to represent this product as a DescriptorSystem with non-unit descriptor matrix if result is inaccurate.") P*tf('s')^3

@test_throws ControlSystemsBase.ImproperException P*tf('s')^5
# Discrete

P = c2d(P, 0.01, :fwdeuler) # to get poles exactly at 1
C = tf('z', 0.01)-1
@test norm(minreal(P*C - C*P)) < 1e-10
mp = minreal(minreal(tf(P)*C) - P*C)
@test hinfnorm(mp)[1] < 1e-10
@inferred P*C

# @test_logs (:warn,"Possible numerical instability detected: Multiplication of a statespace system and a non-proper transfer function may result in numerical inaccuracy. Verify result carefully, and consider making use of DescriptorSystems.jl to represent this product as a DescriptorSystem with non-unit descriptor matrix if result is inaccurate.") P*C^3

@test_throws ControlSystemsBase.ImproperException P*C^5

end
