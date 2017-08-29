using ControlSystems
using Base.Test
import Base.isapprox
include("framework.jl")

my_tests = ["test_statespace",
            "test_transferfunction",
            "test_generalizedtf",
            "test_zpk",
            "test_analysis",
            "test_connections",
            "test_discrete",
            "test_linalg",
            "test_simplification",
            "test_freqresp",
            "test_synthesis",
            "test_matrix_comps",
            "test_timeresp",
            "test_conversion"]

run_tests(my_tests)

# try
#     Pkg.installed("ControlExamplePlots")
#     @testset "test_plots" include("test_plots.jl")
# catch
#     warn("The unregistered package ControlExamplePlots is currently needed to test plots, install using:
#     Pkg.clone(\"https://github.com/JuliaControl/ControlExamplePlots.jl.git\")")
# end
