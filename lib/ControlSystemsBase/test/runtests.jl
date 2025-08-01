using ControlSystemsBase
using Test, LinearAlgebra, Random
import Base.isapprox        # In framework and test_synthesis
import SparseArrays: sparse # In test_matrix_comps
import Polynomials: conv            # In test_conversion and test_synthesis
using Aqua
@testset "Aqua" begin
    Aqua.test_all(ControlSystemsBase;
        ambiguities = false, # causes 100s of hits in all dependencies
        stale_deps = true,  # Aqua complains about itself https://github.com/JuliaTesting/Aqua.jl/issues/78 
        project_toml_formatting = false, # https://github.com/JuliaTesting/Aqua.jl/issues/105
    )
end


include("framework.jl")

my_tests = [
            "test_result_types",
            "test_timeevol",
            "test_statespace",
            "test_staticsystems",
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
            "test_synthesis",
            "test_pid_design",
            "test_partitioned_statespace",
            "test_delayed_systems",
            "test_hammerstein_wiener",
            "test_demo_systems",
            "test_autovec",
            "test_plots",
            "test_dsp",
            "test_implicit_diff",
            "test_rootlocus",
            "test_root_locus_matrix",
            ]

@testset "All Tests" begin
    println("Testing code")
    _t0 = time()
    run_tests(my_tests)
    println("Ran all code tests in $(round(time()-_t0, digits=2)) seconds")
end
