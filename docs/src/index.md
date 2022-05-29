# ControlSystems.jl Manual

```@meta
CurrentModule = ControlSystems
DocTestFilters = [
    r"StateSpace.+?\n"
    r"HeteroStateSpace.+?\n"
    r"TransferFunction.+?\n"
    r"DelayLtiSystem.+?\n"
    r"┌ Warning: Keyword argument hover.+\n*.+\n*" # remove next line as well
    r"\[ Info: Precompiling.+\n*"
]
nyquistplot(ssrand(1,1,1))
```

## Guide

```@contents
Pages = ["man/introduction.md", "man/creating_systems.md", "man/numerical.md", "man/differences.md"]
Depth = 1
```

## Examples
```@contents
Pages = ["examples/example.md", "examples/ilc.md", "examples/smith_predictor.md"]
Depth = 2
```

## Functions

```@contents
Pages = ["lib/constructors.md",  "lib/analysis.md", "lib/syntheis.md", "lib/timefreqresponse.md", "lib/plotting.md"]
Depth = 1
```

## Ecosystem
The JuliaControl and surrounding ecosystem contains a few additional packages that may be of interest
- [RobustAndOptimalControl.jl](https://github.com/JuliaControl/RobustAndOptimalControl.jl)
- [SymbolicControlSystems.jl](https://github.com/JuliaControl/SymbolicControlSystems.jl)
- [ControlSystemIdentification.jl](https://github.com/baggepinnen/ControlSystemIdentification.jl)

See also [the paper](https://portal.research.lu.se/en/publications/controlsystemsjl-a-control-toolbox-in-julia) describing the toolbox

> Bagge Carlson, F., Fält, M., Heimerson, A., & Troeng, O. (2021). ControlSystems.jl: A Control Toolbox in Julia. In 2021 60th IEEE Conference on Decision and Control (CDC) IEEE Press. https://doi.org/10.1109/CDC45484.2021.9683403

and the introductory Youtube video below.

```@raw html
<iframe style="height: 315px; width: 560px" src="https://www.youtube.com/embed/Fdz2Fsm1aTY" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>
```

