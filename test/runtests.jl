include("framework.jl")
using CustomTest

try
    Pkg.installed("ControlExamplePlots")
catch
    warn("ControlExamplePlots needs to be installed to test plots, adding now:")
    Pkg.clone("https://github.com/JuliaControl/ControlExamplePlots.jl.git")
end

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
            "test_conversion",
            "test_plots"]

runtests(my_tests)
