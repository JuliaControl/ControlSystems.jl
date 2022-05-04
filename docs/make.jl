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
    strict=[
        :doctest, 
        :linkcheck, 
        :parse_error, 
        # Other available options are
        # :autodocs_block, :cross_references, :docs_block, :eval_block, :example_block, :footnote, :meta_block, :missing_docs, :setup_block
    ],
    pages=[
        "Home" => "index.md",
        "Guide" => Any[
            "Introduction" => "man/introduction.md",
            "Creating Systems" => "man/creating_systems.md",
        ],
        "Examples" => Any[
            "Design" => "examples/example.md",
            "Smith predictor" => "examples/smith_predictor.md",
            "Iterative Learning Control (ILC)" => "examples/ilc.md",
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
    ]
)

deploydocs(repo = "github.com/JuliaControl/ControlSystems.jl.git")
