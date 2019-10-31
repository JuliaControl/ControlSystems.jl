# Introduction
## Installation

To install this package simply run
```julia
using Pkg; Pkg.add("ControlSystems")
```

## Basic functions
```@meta
DocTestSetup = quote
    using ControlSystems
    P = tf([1],[1,1])
    T = P/(1+P)
    plotsDir = joinpath(dirname(pathof(ControlSystems)), "..", "docs", "build", "plots")
    mkpath(plotsDir)
    save_docs_plot(name) = Plots.savefig(joinpath(plotsDir,name))
    save_docs_plot(p, name) = Plots.savefig(p, joinpath(plotsDir,name))
end
```
Transfer functions can easily be created using the function `tf(num, den, Ts=0)`, where `num` and `den` are vectors representing the numerator and denominator of a rational function. See [`tf`](@ref) or the section "Creating Transfer Functions" for more info. These functions can then be connected and modified using the operators `+,-,*,/` and functions like [`append`](@ref).

Example:
```jldoctest INTRO
P = tf([1.0],[1,1])
T = P/(1+P)

# output

TransferFunction{Continuous, ControlSystems.SisoRational{Float64}}
    1.0s + 1.0
-------------------
1.0s^2 + 3.0s + 2.0

Continuous-time transfer function model
```

Notice that the poles are not canceled automatically, to do this, the function `minreal` is available
```jldoctest INTRO
minreal(T)

# output

TransferFunction{Continuous, ControlSystems.SisoRational{Float64}}
   1.0
----------
1.0s + 2.0

Continuous-time transfer function model
```
## Plotting
Plotting requires some extra care. The ControlSystems package is using `RecipesBase.jl` ([link](https://github.com/JuliaPlots/RecipesBase.jl)) as interface to generate all the plots. This means that it is up to the user to choose a plotting library that supports `RecipesBase.jl`, a suggestion would be `Plots.jl` with which the user is also able to freely choose a back-end. The plots in this manual are generated using `Plots.jl` with the `GR` backend. If you have several back-ends for plotting then you can select the one you want to use with the corresponding `Plots` call (for `GR` this is `Plots.gr()`, some alternatives are `pyplot(), plotly(), pgfplots()`). A simple example where we generate a plot and save it to a file is
```jldoctest; output=false
using Plots

fig = bodeplot(tf(1,[1,2,1]))

savefig(fig, "myfile.svg")

save_docs_plot(fig, "intro_bode.svg") # hide

# output

```

![](../../plots/intro_bode.svg)

More examples of plots are provided in [Plotting](@ref).