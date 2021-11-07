# ControlSystems.jl Manual

```@meta
CurrentModule = ControlSystems
DocTestFilters = [
    r"StateSpace.+?\n"
    r"HeteroStateSpace.+?\n"
    r"TransferFunction.+?\n"
    r"DelayLtiSystem.+?\n"
    r"┌ Warning: Keyword argument hover.+?\n.+?\n" # remove next line as well
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
