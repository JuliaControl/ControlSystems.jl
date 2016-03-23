using Lapidary, ControlSystems, Plots

makedocs(modules=[ControlSystems])

include("src/makeplots.jl")

deploydocs(
    repo = "github.com/JuliaControl/ControlSystems.jl.git",
    latest = "tests",
    julia = "0.4"
)
