using ControlSystems
using Test, LinearAlgebra, Random
using Documenter            # For doctests
import Base.isapprox        # In framework and test_synthesis
import SparseArrays: sparse # In test_matrix_comps
import DSP: conv            # In test_conversion and test_synthesis
include("framework.jl")

# Local definition to simplify tests
eye_(n) = Matrix{Int64}(I, n, n)

my_tests = [
            "test_statespace",
            "test_transferfunction",
            "test_delayed_systems",
            "test_zpk",
            "test_promotion",
            "test_connections",
            "test_discrete",
            "test_conversion",
            "test_complex",
            "test_linalg",
            "test_simplification",
            "test_freqresp",
            "test_timeresp",
            "test_analysis",
            "test_matrix_comps",
            "test_lqg",
            "test_synthesis",
            "test_partitioned_statespace",
            "test_demo_systems",
            ]

@testset "All Tests" begin
    run_tests(my_tests)

    _t0 = time()
    @testset "test_plots_in_docs" begin
        println("Test plots in docs")
        Plots.default(show=false)
        makeplotsfile = joinpath(dirname(pathof(ControlSystems)), "..", "docs", "src", "makeplots.jl")
        include(makeplotsfile)
        # Test that generation of plots needed for documentation is working
        makePlots()
    end
    _t1 = time()
    println("Ran test_plots_in_docs in $(round(_t1-_t0, digits=2)) seconds")


    _t0 = time()
    @testset "Test doctests" begin
        println("test_doctests")
        # Test expected output from examples in documentation
        doctest(ControlSystems)
    end
    _t1 = time()
    println("Ran test_plots_in_docs in $(round(_t1-_t0, digits=2)) seconds")
end
