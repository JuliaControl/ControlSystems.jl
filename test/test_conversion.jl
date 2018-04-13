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

#Test TransferFunction -> StateSpace
#ss is not unique so go back and forth
@test f11 ≈ tf(ss(f12)) atol=1e-15
@test f21 ≈ tf(ss(f22)) atol=1e-15
end
