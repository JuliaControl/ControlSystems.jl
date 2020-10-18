@testset "Test history buffer" begin
    buff = ControlSystems.DDEBuffer(Float64, [1.0,2.0])
    # Getting from empty
    @test_throws ArgumentError buff(0.0, 1)
    @test_throws BoundsError buff(0.0, 3)

    # Push 1
    push!(buff, 1.0, 1.0, 1)
    # Pushing invalid index
    @test_throws BoundsError push!(buff, 2.0, 1.0, 3)
    # Pushing before buffer
    @test_throws AssertionError push!(buff, 0.5, 1.0, 1)

    @test buff(1.0, 1) == 1.0
    # Reading before buffer
    @test_throws AssertionError buff(0.9, 1)
    # Reading empty
    @test_throws ArgumentError buff(0.9, 2)
    # Reading invalid index
    @test_throws BoundsError buff(0.9, 3)

    # Extrapolating
    @test buff(1.5, 1) == 1.0

    push!(buff, 2.0, 2.0, 1)
    # Interpolating
    @test buff(1.5, 1) ≈ 1.5 atol=1e-14

    # Push index 2
    push!(buff, 3.0, 5.0, 2)
    @test buff(3.0, 2) == 5.0
    @test buff(3.1, 2) == 5.0
    push!(buff, 3.1, 6.0, 2)
    @test buff(3.0, 2) == 5.0
    @test buff(3.1, 2) == 6.0
    @test buff(3.09, 2) ≈ 5.9 atol=1e-14

    buff = ControlSystems.DDEBuffer(Float64, [1.0])

    ts = 0:0.1:3

    for t in ts
        push!(buff, t, sin(t), 1)
    end

    xs = Float64[]
    ys = Float64[]

    for i in 1:1000
        t = rand()*3
        y = buff(t,1)
        push!(xs, t)
        push!(ys, y)
    end
    # Roughly maximum error for linear interpolation
    @test maximum(sin.(xs)- ys) < 0.0013
end

@testset "Test history function" begin
    f(t, indx) = indx == 1 ? sin(t) : cos(t)
    hf = ControlSystems.HistoryFunction(Float64, f, 0.0, [1.0, 0.1])

    # Reading before t0  
    @test hf(-1.0, 1) == sin(-1)
    # Reading at t0 and forward not allowd
    @test_throws ArgumentError hf(0.0, 1)

    hf[0.0, 1] = 1.0
    @test hf(0.0, 1) == 1.0
    @test hf(0.1, 1) == 1.0
    hf[0.1, 1] = 2.0
    @test hf(0.1, 1) == 2.0
    @test hf(0.05, 1) ≈ 1.5 atol=1e-14

    # Inserting at old time
    @test_throws AssertionError hf[0.05,1] = 1.0

    # No data at index 2
    @test_throws ArgumentError hf(0.05, 2)

    for t in 0.2:0.01:0.99
        hf[t,1] = t
        hf[t,2] = 2*t
    end

    @test hf(0.21,1) == 0.21
    # Check old value
    @test hf(0.0, 1) == 1
    @test hf(0.99, 1) == 0.99
    # Test extrapolation
    @test hf(1.01, 1) == 0.99

    # Test old value for 2
    @test hf(0.89, 2) == 0.89*2
    # This should be dropped
    @test_throws AssertionError hf(0.87, 2)
    # This should just barely be in
    @test hf(0.88, 2) == 0.88*2
    # Add within epsilon away from being able to read 0.99-tau
    hf[0.99+1e-14,2] = 1.0
    # We should still be able to read
    @test hf(0.88, 2) == 0.88*2
    # Some interpolation
    @test hf(0.985,1) ≈ 0.985 atol=1e-14
    @test hf(0.985,2) ≈ 0.985*2 atol=1e-14
end