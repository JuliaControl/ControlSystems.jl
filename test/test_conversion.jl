@testset "test_conversion" begin

f11 = tf(ss([-1 0;1 1],[1;0],[1 1],0))
f12 = tf([1,0],[1,0,-1])

#Direct term was previously missing
f21 = tf(ss([-1 0;1 1],[1;0],[1 1],1))
f22 = tf([1,1,-1],[1,0,-1])

f31 = tf(ss([0 1;1 0],[0;1],[1 1],0))
f32 = tf([1,1],[1,0,-1])

#Test ss2tf
@test f11 ≈ f12 atol=1e-15
@test f21 ≈ f22 atol=1e-15
@test f31 ≈ f32 atol=1e-15

#Test complex
A = [1.0 + im 1; 0 -2-3im]
B = [0;2]
C = [1 0]
D = zeros(1,1);
sys = ss(A,B,C,D)
sys2 = convert(TransferFunction, sys)
w = 10.0 .^ range(-2,2,50)
@test freqresp(sys, w) ≈ freqresp(sys2, w)
csort = v -> sort(v, lt = (x,y) -> abs(x) < abs(y))
@test csort(pole(zpk(sys2))) ≈ [1+im, -2.0-3im]

#Test TransferFunction -> StateSpace
#ss is not unique so go back and forth
@test f11 ≈ tf(ss(f12)) atol=1e-15
@test f21 ≈ tf(ss(f22)) atol=1e-15

## Test some rather trivial conversions of numeric types

b = 1.5
D22 = [1.0 2.0; 3.0 4.0]

@test convert(StateSpace{Float64,Matrix{Float64}}, D22) == ss(D22)
@test convert(StateSpace{Float64,Matrix{Float64}}, b) == ss(b)

@test convert(TransferFunction{ControlSystems.SisoRational{Float64}}, D22) == tf(D22)
@test convert(TransferFunction{ControlSystems.SisoRational{Float64}}, b) == tf(b)

@test convert(TransferFunction{ControlSystems.SisoZpk{Float64,Complex128}}, D22) == zpk(D22)
@test convert(TransferFunction{ControlSystems.SisoZpk{Float64,Complex128}}, b) == zpk(b)


# Error, not proper
@test_throws ErrorException ss(tf([1,0],[1]))

end
