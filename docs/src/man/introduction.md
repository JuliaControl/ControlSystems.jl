# Introduction
## Installation

To install the full set of features, simply run
```julia
using Pkg; Pkg.add("ControlSystems")
```
For workflows that do not require continuous-time simulation, you may instead opt to install the *much* lighter package ControlSystemsBase.jl
```julia
using Pkg; Pkg.add("ControlSystemsBase")
```
ControlSystemsBase contains all functionality of ControlSystems except continuous-time simulation and root locus, and is *considerably* faster to load and precompile. To enjoy the faster pre-compilation, do not even install ControlSystems since this will cause pre-compilation of OrdinaryDiffEq, which can take several minutes.

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
State-space systems can be created using the function [`ss`](@ref) and transfer functions can be created using the function `tf(num, den)` or `tf(num, den, Ts)`, where `num` and `den` are vectors representing the numerator and denominator of a rational function and `Ts` is the sample time for a discrete-time system. See [`tf`](@ref) or the section [Creating Systems] for more info. These functions can then be connected and modified using the operators `+,-,*,/` and functions like [`append`](@ref).

Example:
```jldoctest INTRO
P = tf([1.0],[1,1])
T = P/(1+P)

# output

TransferFunction{Continuous, ControlSystemsBase.SisoRational{Float64}}
    1.0s + 1.0
-------------------
1.0s^2 + 3.0s + 2.0

Continuous-time transfer function model
```

Notice that the poles are not canceled automatically, to do this, the function `minreal` is available
```jldoctest INTRO
minreal(T)

# output

TransferFunction{Continuous, ControlSystemsBase.SisoRational{Float64}}
   1.0
----------
1.0s + 2.0

Continuous-time transfer function model
```
or use `feedback(P)` to get a minimal realization directly (recommended):
```@example INTRO
using ControlSystems # hide
P = tf([1.0],[1,1]) # hide
feedback(P) # Equivalent to P/(1+P)
```

!!! note "Numerical accuracy"
    Transfer functions represent systems using polynomials and may have poor numerical properties for high-order systems. Well-balanced state-space representations are often better behaved. See [Performance considerations](@ref) for more details.

## Plotting
The ControlSystems package is using `RecipesBase.jl` ([link](https://github.com/JuliaPlots/RecipesBase.jl)) as interface to generate all the plots. This means that it is up to the user to choose a plotting library that supports `RecipesBase.jl`, a suggestion would be `Plots.jl` with which the user is also able to freely choose a back-end. The plots in this manual are generated using `Plots.jl` with the `GR` backend. If you have several back-ends for plotting then you can select the one you want to use with the corresponding `Plots` call (for `GR` this is `Plots.gr()`, some alternatives are `pyplot(), plotly(), pgfplots()`). A simple example where we generate a plot and save it to a file is shown below.

More examples of plots are provided in [Plotting functions](@ref).

```@example INTRO
using Plots
bodeplot(tf(1,[1,2,1]))
```

