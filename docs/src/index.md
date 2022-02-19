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
Pages = ["man/introduction.md", "man/creating_systems.md"]
Depth = 1
```

## Examples
```@contents
Pages = ["examples/example.md"]
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
