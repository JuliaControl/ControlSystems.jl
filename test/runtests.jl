using ControlSystems
using Test, LinearAlgebra, Random
import Base.isapprox        # In framework and test_synthesis
import SparseArrays: sparse # In test_matrix_comps
import DSP: conv            # In test_conversion and test_synthesis
include("framework.jl")

# Local definition to make sure we get warnings if we use eye
eye_(n) = Matrix{Int64}(I, n, n)

my_tests = [
            "test_timeevol",
            "test_statespace",
            "test_transferfunction",            
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
            "test_pid_design",
            "test_partitioned_statespace",
            "test_delayed_systems",
            "test_demo_systems",
            "test_autovec",
            "test_plots"
            ]

@testset "All Tests" begin
    println("Testing code")
    _t0 = time()
    run_tests(my_tests)
    println("Ran all code tests in $(round(time()-_t0, digits=2)) seconds")
end
