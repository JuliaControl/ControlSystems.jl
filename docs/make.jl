# Set plot globals
ENV["PLOTS_TEST"] = "true"
ENV["GKSwstype"] = "100"

using Documenter, ControlSystems, Plots, LinearAlgebra, DSP
import GR # Bug with world age in Plots.jl, see https://github.com/JuliaPlots/Plots.jl/issues/1047

dir = joinpath(dirname(pathof(ControlSystems)), "..")
cd(dir)
include(joinpath(dir, "docs", "src", "makeplots.jl"))

println("Making plots for docs")
makePlots()

println("Making docs")
makedocs(modules=[ControlSystems],
    format=Documenter.HTML(),
    sitename="ControlSystems.jl",
    pages=[
        "Home" => "index.md",
        "Examples" => Any[
            "Design" => "examples/example.md",
        ],
        "Guide" => Any[
            "Introduction" => "man/introduction.md",
            "Creating Systems" => "man/creating_systems.md",
        ],
        "Functions" => Any[
            "Constructors" => "lib/constructors.md",
            "Analysis" => "lib/analysis.md",
            "Synthesis" => "lib/synthesis.md",
            "Time and Frequency response" => "lib/timefreqresponse.md",
            "Plotting" => "lib/plotting.md",
        ],
    ]
    )

deploydocs(repo = "github.com/JuliaControl/ControlSystems.jl.git")
