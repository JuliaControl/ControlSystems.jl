using Control
include("systems.jl")
include("framework.jl")

my_tests = ["test_statespace",
            "test_linalg"]

runtests(my_tests)
