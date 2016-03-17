include("framework.jl")
using CustomTest

my_tests = ["test_statespace",
            "test_transferfunction",
            "test_zpk",
            "test_analysis",
            "test_connections",
            "test_discrete",
            "test_linalg",
            "test_simplification",
            "test_freqresp"]

runtests(my_tests)
