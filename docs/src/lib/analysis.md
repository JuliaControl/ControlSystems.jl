```@index
Pages = ["analysis.md"]
```

For robust analysis, see [RobustAndOptimalControl.jl](https://juliacontrol.github.io/RobustAndOptimalControl.jl/dev/#System-analysis).

# Analysis

```@autodocs
Modules = [ControlSystems, ControlSystemsBase]
Pages   = [
         libpath*"/analysis.jl", 
         libpath*"/matrix_comps.jl", 
         libpath*"/conversion.jl"
        ]
Order   = [:function, :type]
Private = false
```

## Videos

Basic usage of robustness analysis with JuliaControl
```@raw html
<iframe style="height: 315px; width: 560px" src="https://www.youtube.com/embed/zTW4mlWNumo" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>
```