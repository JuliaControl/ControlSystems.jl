using ControlSystems
using Test, LinearAlgebra, Random
import Base.isapprox
import SparseArrays: sparse
import DSP: conv
include("framework.jl")

# Local definition to make sure we get warings if we use eye
eye_(n) = Matrix{Int64}(I, n, n)

my_tests = ["test_statespace",
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
            "test_synthesis"]


# try
#     Pkg.installed("ControlExamplePlots")
#     push!(my_tests, "test_plots")
# catch
#     warn("The unregistered package ControlExamplePlots is currently needed to test plots, install using:
#     Pkg.clone(\"https://github.com/JuliaControl/ControlExamplePlots.jl.git\")")
# end

run_tests(my_tests)
