using Test, ControlSystems.ReferenceGenerators 

# PWC tests
@testset "PWCGenerator Tests" begin 
    pwc1 = PWCGenerator([1, 2], [1, 0])
    pwc2 = PWCGenerator([-1, 0, 1], [1, -1, 0]) 

    @test pwc1(1.5) == 1
    @test pwc1(0.5) == 0
    @test pwc1(2.5) == 0

    @test pwc1.([0.99, 1.00, 1.01]) == [0, 1, 1]
    @test pwc1.([1.99, 2.00, 2.01]) == [1, 0, 0]

    @test pwc2.([-1.01, -1.00, -0.99]) == [0, 1, 1]
    @test pwc2.([-0.01, 0., 0.01])     == [1, -1, -1]
    @test pwc2.([0.99, 1.00, 1.01])    == [-1, 0, 0]
    @test pwc2.([-10, -0.5, 0.5, 10])  == [0, 1, -1, 0]
end

# Chirp tests
@testset "ChirpGenerator Tests" begin
    c1       = ChirpGenerator()
    test1(t) = sin(t^2)
    c2       = ChirpGenerator(t -> t^2 + 10)
    test2(t) = sin((t^2+10)*t)
    c3       = ChirpGenerator(t -> t + 10)

    @test c1(0) == 0
    @test c1(1) == sin(1)

    t = 0:0.001:10
    @test c1.(t) == test1.(t)
    @test c2.(t) == test2.(t)
end 

# Repeating Sequence Tests
@testset "RSGenerator Tests" begin
    rs1      = RSGenerator(t -> sin(t*(2pi)), 1)
    test1(t) = sin(t*(2pi))
    pwc1     = PWCGenerator([1, 2], [1, 0])
    rs2      = RSGenerator(pwc1, 0.5)
    rs3      = RSGenerator(t -> sin(t*(2pi)), 2)

    @test rs2(0.5) == 0
    @test rs2(1.5) == 1
    @test rs2(2.5) == 0
    @test rs2(3.5) == 1

    t = 0:0.01:4
    @test norm(rs1.(t) - test1.(t)) < 1e-13
    @test rs2.(0:0.01:2) == rs2.(2:0.01:4)
    @test norm(rs3.(0:0.01:0.5) - rs3.(1.5:0.01:2)) < 1e-13
end 

# Sawtooth Tests
@testset "SawtoothGenerator Tests" begin
    s1 = SawtoothGenerator(1)
    s2 = SawtoothGenerator(1., 0.5)
    s3 = SawtoothGenerator(2)

    @test s2(0) == 0.5

    @test norm(s1.(0:0.01:1) - s1.(2:0.01:3)) < 1e-13
    @test norm(s1.(-1:0.01:0) - s1.(2:0.01:3)) < 1e-13
    @test norm(s2.(0:0.01:1) - s2.(2:0.01:3)) < 1e-13
    @test norm(s3.(0:0.01:0.5) - s3.(1:0.01:1.5)) < 1e-13
end 

# Square Tests
@testset "SquareGenerator Tests" begin
    s1 = SquareGenerator(1)
    s2 = SquareGenerator(1, 5)
    s3 = SquareGenerator(4, 5)

    @test s2(0.00) == 5
    @test s2(0.25) == 5
    @test s2(0.50) == -5
    @test s2(0.75) == -5
    @test s2(1.00) == 5

    @test s1.A == 1
    @test norm(s1.(0:0.01:1) - s1.(2:0.01:3)) < 1e-13
    @test s2.([-0.1, 0, 0.2, 0.5, 0.7, 1]) == [-5, 5, 5, -5, -5, 5]
    @test norm(s1.(0:0.01:0.25) - s1.(2:0.01:2.25)) < 1e-13
end 

# Gaussian Noise Tests
@testset "GNGenerator Tests" begin
    mean(x) = sum(x) / length(x)
    gn1 = GNGenerator(1)
    gn2 = GNGenerator(1, 1)

    t = zeros(1_000_000)
    @test abs(mean(gn1.(t))) < 1e-2
    @test abs(mean(gn2.(t)) - 1)< 1e-2
end 

# Noise Tests
@testset "NoiseGenerator Tests" begin
    using Distributions
    n1 = NoiseGenerator(Poisson(1))
    #TODO: Add more distributions?

    t = zeros(1_000_000)
    @test abs(mean(n1.(t)) - 1)< 1e-2
end 

# RG Collector Tests
@testset "RGCollector Tests" begin
    s1    = SawtoothGenerator(1)
    pwc1  = PWCGenerator([1, 2], [1, 0])
    f1(t) = sin(4pi*t)
    rgc   = RGCollector(s1, pwc1, f1)

    t     = 0:0.01:4
    vals  = rgc.(t)

    @test rgc(0) == [0, 0, 0]

    #@test norm(vals[1] - s1.(t)) < 1e-13
    #@test norm(vals[2] - pwc1.(t)) < 1e-13
    #@test norm(vals[3] - f1.(t)) < 1e-13
end 
