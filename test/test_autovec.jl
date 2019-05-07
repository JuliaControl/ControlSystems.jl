@testset "test_autovec" begin

# Test creating simple function and check output
ControlSystems.@autovec (1, 3) 4 f() = (ones(3, 1, 1), ones(3), ones(1, 3), ones(3))
@test f() == (ones(3, 1, 1), ones(3), ones(1, 3), ones(3))
@test fv() == (ones(3), ones(3), ones(3), ones(3))

# Check output of bodev and make sure dimensions are correct
sys = tf([1], [1, 2])
mag, phase, w = bodev(sys)
@test size(mag) == size(phase) == size(w)

end
