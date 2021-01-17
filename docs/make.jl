# Set plot globals
ENV["PLOTS_TEST"] = "true"
ENV["GKSwstype"] = "nul"

using Documenter, ControlSystems, Plots, LinearAlgebra, DSP
#import GR # Bug with world age in Plots.jl, see https://github.com/JuliaPlots/Plots.jl/issues/1047
gr()
default(show=false, size=(800,450))

dir = joinpath(dirname(pathof(ControlSystems)), "..")
cd(dir)

# Copied from Documenter/src/Document.jl, modified to remove # hide lines
Markdown = Documenter.Documents.Markdown
function Documenter.Documents.doctest_replace!(block::Markdown.Code)
    startswith(block.language, "jldoctest") || return false
    # suppress output for `#output`-style doctests with `output=false` kwarg
    if occursin(r"^# output$"m, block.code) && occursin(r";.*output\h*=\h*false", block.language)
        input = first(split(block.code, "# output\n", limit = 2))
        block.code = rstrip(input)
    end
    # Remove # hide lines
    block.code = Documenter.Expanders.droplines(block.code)
    # correct the language field
    block.language = occursin(r"^julia> "m, block.code) ? "julia-repl" : "julia"
    return false
end

println("Making docs")
makedocs(modules=[ControlSystems],
    format=Documenter.HTML(),
    sitename="ControlSystems.jl",
    #strict=true,
    pages=[
        "Home" => "index.md",
        "Guide" => Any[
            "Introduction" => "man/introduction.md",
            "Creating Systems" => "man/creating_systems.md",
        ],
        "Examples" => Any[
            "Design" => "examples/example.md",
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
