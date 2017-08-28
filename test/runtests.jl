include("framework.jl")
using CustomTest
using Base.Test

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

try
    Pkg.installed("ControlExamplePlots")
    push!(my_tests, "test_plots")
catch
    warn("The unregistered package ControlExamplePlots is currently needed to test plots, install using:
            Pkg.clone(\"https://github.com/JuliaControl/ControlExamplePlots.jl.git\")")
end

run_tests(my_tests)
