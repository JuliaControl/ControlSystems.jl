using ControlSystems
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

end