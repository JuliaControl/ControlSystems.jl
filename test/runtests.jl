using ControlSystems
using Test, LinearAlgebra, Random
using Documenter # For doctests
import Base.isapprox
import SparseArrays: sparse
import DSP: conv
include("framework.jl")

# Local definition to make sure we get warings if we use eye
eye_(n) = Matrix{Int64}(I, n, n)

my_tests = [
            "test_statespace",
            # "test_transferfunction",
            # "test_delayed_systems",
            # "test_zpk",
            # "test_promotion",
            # "test_connections",
            # "test_discrete",
            # "test_conversion",
            # "test_complex",
            # "test_linalg",
            # "test_simplification",
            # "test_freqresp",
            # "test_timeresp",
            # "test_analysis",
            # "test_matrix_comps",
            # "test_lqg",
            # "test_synthesis",
            # "test_partitioned_statespace",
            # "test_demo_systems",
            ]

@testset "All Tests" begin
    run_tests(my_tests)
    @testset "Test Doc plots" begin
        Plots.default(show=false)
        include(joinpath(dirname(pathof(ControlSystems)), "..", "docs", "src", "makeplots.jl"))
        makePlots()
    end
    @testset "Doctests" begin
        doctest(ControlSystems, fix=true)
    end
end
