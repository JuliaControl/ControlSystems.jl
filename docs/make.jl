using Documenter, ControlSystems, Plots, LinearAlgebra, DSP
import GR # Bug with world age in Plots.jl, see https://github.com/JuliaPlots/Plots.jl/issues/1047

include("src/makeplots.jl")

makedocs(modules=[ControlSystems], format=:html, sitename="ControlSystems")

# If not running travis, generate the plots here, even if we are not deploying
if get(ENV, "TRAVIS", "") == ""
    makePlots()
end

# Only build plots in travis if we are deploying
# And dont install the dependencies unless we are deploying
function myDeps()
    if get(ENV, "TRAVIS", "") != ""
        println("Installing deploy dependencies")
        makePlots()
    end
end

deploydocs(
    repo = "github.com/JuliaControl/ControlSystems.jl.git",
    deps = myDeps
)
