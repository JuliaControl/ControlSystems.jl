include("framework.jl")
using CustomTest

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

runtests(my_tests)
