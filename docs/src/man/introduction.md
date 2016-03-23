# Introduction

## Installation

To install this package simply run
## Installation
```julia
Pkg.add("ControlSystems")
```

## Basic functions

## Plotting
Plotting requires some extra care. The ControlSystems package is using `Plots.jl` as interface to generate all the plots. This means that the user is able to freely choose back-end. The plots in this manual are generated using `PyPlot`. If you have several back-ends for plotting then you can select the one you want to use with the corresponding `Plots` call (for `PyPlot` this is `Plots.pyplot()`, some alternatives are `immerse(), gadfly(), plotly()`). A simple example where we generate a plot using `immerse` and save it to a file is
```julia
using ControlSystems
Plots.immerse()

fig = bodeplot(tf(1,[1,2,1]))

Plots.savefig(fig, "myfile.svg")
```