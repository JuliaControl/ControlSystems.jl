# Set plot globals
ENV["PLOTS_TEST"] = "true"
ENV["GKSwstype"] = "nul"

using Documenter, ControlSystems, ControlSystemsBase, Plots, LinearAlgebra
import DSP
# ENV["JULIA_DEBUG"]=Documenter # Enable this for debugging
#import GR # Bug with world age in Plots.jl, see https://github.com/JuliaPlots/Plots.jl/issues/1047
gr()
default(show=false, size=(800,450))

dir = joinpath(dirname(pathof(ControlSystems)), "..")
cd(dir)

const libpath = haskey(ENV, "CI") ? dirname(pathof(ControlSystemsBase)) : "lib/ControlSystemsBase/src"
dirname(pathof(ControlSystemsBase))

# Update doctest outputs with doctest("/home/fredrikb/.julia/dev/ControlSystems/docs", [ControlSystems], fix=true)

println("Making docs")
makedocs(modules=[ControlSystems, ControlSystemsBase],
    format = Documenter.HTML(prettyurls = haskey(ENV, "CI")),
    sitename="ControlSystems.jl",
    pagesonly = true,
    draft = false,
    pages=[
        "Home" => "index.md",
        "Introductory guide" => Any[
            "Introduction" => "man/introduction.md",
            "Creating Systems" => "man/creating_systems.md",
            "Performance considerations" => "man/numerical.md",
            "Noteworthy differences from other languages" => "man/differences.md",
        ],
        "Examples" => Any[
            "Design" => "examples/example.md",
            "Analysis" => "examples/analysis.md",
            "Smith predictor" => "examples/smith_predictor.md",
            "Iterative Learning Control (ILC)" => "examples/ilc.md",
            "Properties of delay systems" => "examples/delay_systems.md",
            "Automatic differentiation" => "examples/automatic_differentiation.md",
            "Tune a controller using experimental data" => "examples/tuning_from_data.md",
            "Analysis of hybrid continuous/discrete systems" => "examples/zoh.md",
        ],
        "Functions" => Any[
            "Constructors" => "lib/constructors.md",
            "Analysis" => "lib/analysis.md",
            "Synthesis" => "lib/synthesis.md",
            "Time and Frequency response" => "lib/timefreqresponse.md",
            "Plotting" => "lib/plotting.md",
            "Nonlinear" => "lib/nonlinear.md",
        ],
        "API" => "api.md",
    ],
    warnonly = [:missing_docs, :cross_references],
)

deploydocs(repo = "github.com/JuliaControl/ControlSystems.jl.git")
