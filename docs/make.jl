using Lapidary, ControlSystems, Plots

makedocs(modules=[ControlSystems])

include("src/makeplots.jl")

deploydocs(
    repo = "github.com/JuliaControl/ControlSystems.jl.git",
    latest = "master",
    julia = "0.4",
    deps = myDeps()
)

function myDeps()
    if get(ENV, "TRAVIS", "") != ""
        run(`pip install --user pygments mkdocs`)
        ENV["PYTHON"]=""
        Pkg.add("PyPlot")
        Pkg.build("PyCall")
    end
end
