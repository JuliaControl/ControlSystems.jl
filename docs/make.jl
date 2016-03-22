using Lapidary, ControlSystems

makedocs()

deploydocs(
    repo = "github.com/JuliaControl/ControlSystems.jl.git",
    latest = "tests",
    julia = "0.4"
)
