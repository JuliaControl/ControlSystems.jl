@testset "test_autovec" begin

# Check output of bode and make sure dimensions are correct
sys = tf([1], [1, 2])
mag, phase, w = bode(sys)
@test size(mag) == size(phase) == (size(w)[1], 1, 1)

# Check output of bodev and make sure dimensions are correct
mag, phase, w = bodev(sys)
@test size(mag) == size(phase) == size(w)

end
