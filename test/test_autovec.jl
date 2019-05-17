module TestAutovec
import ControlSystems: @autovec

@autovec (1, 3) 4 function test1(a::T, b::Int; c=1) where {T<:TransferFunction}
    return ones(3, 1), ones(3), ones(1, 3), ones(3)
end

end

@testset "test_autovec" begin

using Main.TestAutovec

A = [-1.0 -2.0; 0.0 -1.0]
B = [0.0; -2.0]
C = [1.0 0; 0 1.0]
D = [1.0; 0]

sys_simo = ss(A, B, C, D)

@test test1(sys_simo, 3)
@test_throws ArgumentError test1v(sys_simo, 3)


sys = tf([1], [1, 2])

# Check output of bode and make sure dimensions are correct
mag, phase, w = bode(sys)
@test size(mag) == size(phase) == (size(w,1), 1, 1)

# Check output of bodev and make sure dimensions are correct
mag, phase, w = bodev(sys)
@test size(mag) == size(phase) == size(w)

end
