using Lapidary, ControlSystems, Plots
include("src/makeplots.jl")

makedocs(modules=[ControlSystems])

# If not running travis, generate the plots here, even if we are not deploying
if get(ENV, "TRAVIS", "") == ""
    makePlots()
end

# Only build plots in travis if we are deploying
# And dont install the dependencies unless we are deploying
function myDeps()
    if get(ENV, "TRAVIS", "") != ""
        run(`pip install --user pygments mkdocs`)
        ENV["PYTHON"]=""
        Pkg.add("PyPlot")
        Pkg.build("PyCall")
        makePlots()
    end
end

deploydocs(
    repo = "github.com/JuliaControl/ControlSystems.jl.git",
    latest = "master",
    julia = "0.4",
    deps = myDeps()
)
