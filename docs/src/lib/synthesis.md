```@index
Pages = ["synthesis.md"]
```

# Synthesis

For ``H_\infty`` and ``H_2`` synthesis as well as more advanced LQG design, see [RobustAndOptimalControl](https://juliacontrol.github.io/RobustAndOptimalControl.jl/dev/#H_\\infty-and-H_2-design).

```@autodocs
Modules = [ControlSystems, ControlSystemsBase]
Pages   = [
        libpath*"/synthesis.jl",
        libpath*"/discrete.jl",
        libpath*"/types/lqg.jl",
        libpath*"/pid_design.jl",
        libpath*"/simplification.jl",
        libpath*"/connections.jl",
        libpath*"/sensitivity_functions.jl",
        libpath*"/utilities.jl"
    ]
Order   = [:function, :type]
Private = false
```
