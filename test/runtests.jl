using ControlSystems
using Test, LinearAlgebra, Random
using Aqua
# @testset "Aqua" begin
#     Aqua.test_all(ControlSystems;
#         ambiguities = false, # casues 100s of hits in all dependencies
#         stale_deps = true,  # Aqua complains about itself https://github.com/JuliaTesting/Aqua.jl/issues/78 
#     )
# end

@testset "ControlSystems" begin
    @info "Testing ControlSystems"

    @testset "timeresp" begin
        @info "Testing timeresp"
        include("test_timeresp.jl")
    end

    @testset "delay_timeresp" begin
        @info "Testing delay_timeresp"
        include("test_delay_timeresp.jl")
    end

    @testset "nonlinear_timeresp" begin
        @info "Testing nonlinear_timeresp"
        include("test_nonlinear_timeresp.jl")
    end

    @testset "rootlocus" begin
        @info "Testing rootlocus"
        include("test_rootlocus.jl")
    end
end