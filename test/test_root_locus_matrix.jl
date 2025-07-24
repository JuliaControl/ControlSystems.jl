using ControlSystems
using ControlSystems: getpoles
using Test
using Plots

@testset "Root Locus with Matrix Feedback" begin
    
    @testset "State Feedback Matrix" begin
        # Test with double mass model
        P = DemoSystems.double_mass_model()
        
        # Design a simple state feedback matrix
        K_state = [1.0 0.5 0.2 0.1]  # nu x nx matrix (1x4)
        
        # Test that getpoles works with matrix feedback
        roots, K_values = getpoles(P, K_state)
        
        # Basic sanity checks
        @test size(roots, 2) == size(P.A, 1)  # Should have nx poles
        @test length(K_values) == size(roots, 1)  # K_values should match number of steps
        @test eltype(K_values) <: AbstractMatrix  # K_values should be matrices
        
        # Test that initial K_value is zero matrix
        @test all(K_values[1] .≈ 0.0)
        
        # Test that final K_value is the input matrix
        @test K_values[end] ≈ K_state
        
        # Test poles are continuous (no sudden jumps)
        # The poles should change smoothly
        pole_diffs = diff(roots, dims=1)
        max_diff = maximum(abs, pole_diffs)
        @test max_diff < 1.0  # Poles shouldn't jump too much between steps
        
        # Test that rlocus works with matrix K
        result = rlocus(P, K_state)
        @test result isa ControlSystemsBase.RootLocusResult
        @test size(result.roots) == size(roots)
        
        # Test plotting (should not error)
        plt = plot(result)
        @test plt isa Plots.Plot
    end
    
    @testset "Output Feedback Matrix" begin
        P = DemoSystems.double_mass_model(outputs=1:2)
        
        # Design output feedback matrix
        # P has 2 outputs, 1 input, so K should be 1x2
        K_output = [1.0 0.5]  # nu x ny matrix (1x2)
        
        # Test output feedback
        roots, K_values = getpoles(P, K_output; output=true)
        
        # Basic checks
        @test size(roots, 2) == size(P.A, 1)  # Should have nx poles
        @test length(K_values) == size(roots, 1)
        @test eltype(K_values) <: AbstractMatrix
        
        # Test rlocus with output feedback
        result = rlocus(P, K_output; output=true)
        @test result isa ControlSystemsBase.RootLocusResult
        
        # Test plotting
        plt = plot(result)
        @test plt isa Plots.Plot
    end
    
    @testset "Error Handling" begin
        P = DemoSystems.double_mass_model()
        
        # Test wrong dimensions for state feedback
        K_wrong = [1.0 0.5 0.2]  # Wrong size (1x3 instead of 1x4)
        @test_throws ErrorException getpoles(P, K_wrong)
        
        # Test wrong dimensions for output feedback
        K_wrong_output = [1.0 0.5 0.2]  # Wrong size (1x3 instead of 1x2)
        @test_throws ErrorException getpoles(P, K_wrong_output; output=true)
        
        # Test wrong number of rows
        K_wrong_rows = [1.0 0.5 0.2 0.1; 2.0 1.0 0.4 0.2]  # 2x4 instead of 1x4
        @test_throws ErrorException getpoles(P, K_wrong_rows)
    end
    

    @testset "Adaptive Step Size" begin
        P = DemoSystems.double_mass_model()
        K_state = [1.0 0.5 0.2 0.1]
        
        # Test with different tolerances
        roots_loose, K_loose = getpoles(P, K_state; tol=1e-1)
        roots_tight, K_tight = getpoles(P, K_state; tol=1e-4)
        
        # Tighter tolerance should result in more steps
        @test length(K_tight) >= length(K_loose)
        
        # Both should start and end at the same points
        @test roots_loose[1, :] ≈ roots_tight[1, :]  # Same initial poles
        @test roots_loose[end, :] ≈ roots_tight[end, :] atol=1e-3  # Similar final poles
    end
    
end