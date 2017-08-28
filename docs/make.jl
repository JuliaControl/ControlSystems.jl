using Documenter, ControlSystems, Plots
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
        println("Installing deploy dependencies")
        run(`pip install --user pygments mkdocs`)
        makePlots()
    end
end

deploydocs(
    repo = "github.com/JuliaControl/ControlSystems.jl.git",
    latest = "master",
    julia = "0.5",
    deps = myDeps
)
