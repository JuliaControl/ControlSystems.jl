include("framework.jl")
using CustomTest

my_tests = ["test_statespace",
            "test_transferfunction",
            "test_analysis",
            "test_connections",
            "test_discrete",
            "test_linalg"]

runtests(my_tests)
