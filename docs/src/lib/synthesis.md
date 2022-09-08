```@index
Pages = ["synthesis.md"]
```

# Synthesis

For ``H_\infty`` and ``H_2`` synthesis as well as more advanced LQG design, see [RobustAndOptimalControl](https://juliacontrol.github.io/RobustAndOptimalControl.jl/dev/#H_\\infty-and-H_2-design).

```@autodocs
Modules = [ControlSystems, ControlSystemsBase]
Pages   = ["lib/ControlSystemsBase/src/synthesis.jl", "lib/ControlSystemsBase/src/discrete.jl", "lib/ControlSystemsBase/src/types/lqg.jl", "lib/ControlSystemsBase/src/pid_design.jl", "lib/ControlSystemsBase/src/simplification.jl", "lib/ControlSystemsBase/src/connections.jl", "lib/ControlSystemsBase/src/sensitivity_functions.jl", "lib/ControlSystemsBase/src/utilities.jl"]
Order   = [:function, :type]
Private = false
```
