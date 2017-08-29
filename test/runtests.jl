using ControlSystems
using Base.Test
import Base.isapprox
include("framework.jl")

@testset "test_statespace" begin include("test_statespace.jl") end
@testset "test_transferfunction" begin include("test_transferfunction.jl") end
@testset "test_generalizedtf" begin include("test_generalizedtf.jl") end
@testset "test_zpk" begin include("test_zpk.jl") end
@testset "test_analysis" begin include("test_analysis.jl") end
@testset "test_connections" begin include("test_connections.jl") end
@testset "test_discrete" begin include("test_discrete.jl") end
@testset "test_linalg" begin include("test_linalg.jl") end
@testset "test_simplification" begin include("test_simplification.jl") end
@testset "test_freqresp" begin include("test_freqresp.jl") end
@testset "test_synthesis" begin include("test_synthesis.jl") end
@testset "test_matrix_comps" begin include("test_matrix_comps.jl") end
@testset "test_timeresp" begin include("test_timeresp.jl") end
@testset "test_conversion" begin include("test_conversion.jl") end


# try
#     Pkg.installed("ControlExamplePlots")
#     @testset "test_plots" include("test_plots.jl")
# catch
#     warn("The unregistered package ControlExamplePlots is currently needed to test plots, install using:
#     Pkg.clone(\"https://github.com/JuliaControl/ControlExamplePlots.jl.git\")")
# end
