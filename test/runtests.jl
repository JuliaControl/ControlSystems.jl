using ControlSystems
using Test, LinearAlgebra, Pkg
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

if get(ENV, "TRAVIS", "") != ""
    # TODO how to do this without adding in REQUIRE?
    Pkg.add(PackageSpec(url="https://github.com/JuliaControl/ControlExamplePlots.jl", version="0.2.0"))
end

try
    using ControlExamplePlots
    push!(my_tests, "test_plots")
catch
    @warn "The unregistered package ControlExamplePlots is currently needed to test plots, install using:
          Pkg.add(\"https://github.com/JuliaControl/ControlExamplePlots.jl.git\")"
end

run_tests(my_tests)
