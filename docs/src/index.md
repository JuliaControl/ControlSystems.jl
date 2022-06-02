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

### JuliaControl

The JuliaControl and surrounding ecosystem contains a few additional packages that may be of interest
- [RobustAndOptimalControl.jl](https://github.com/JuliaControl/RobustAndOptimalControl.jl) contains more advanced features for LQG design, robust analysis and synthesis, uncertainty modeling, named systems and an interface to [DescriptorSystems.jl](https://github.com/andreasvarga/DescriptorSystems.jl).
- [SymbolicControlSystems.jl](https://github.com/JuliaControl/SymbolicControlSystems.jl) contains basic C-code generation for linear systems.
- [ControlSystemIdentification.jl](https://github.com/baggepinnen/ControlSystemIdentification.jl) is a system-identification toolbox for identification of LTI systems using either time or frequency-domain data.

See also [the paper](https://portal.research.lu.se/en/publications/controlsystemsjl-a-control-toolbox-in-julia) describing the toolbox

> Bagge Carlson, F., Fält, M., Heimerson, A., & Troeng, O. (2021). ControlSystems.jl: A Control Toolbox in Julia. In 2021 60th IEEE Conference on Decision and Control (CDC) IEEE Press. https://doi.org/10.1109/CDC45484.2021.9683403

and the introductory Youtube video below.

```@raw html
<iframe style="height: 315px; width: 560px" src="https://www.youtube.com/embed/Fdz2Fsm1aTY" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>
```

### The wider Julia ecosystem for control
The following is a list of packages from the wider Julia ecosystem that may be of interest.

- [DescriptorSystems.jl](https://github.com/andreasvarga/DescriptorSystems.jl) contains types that represent statespace systems on descriptor form, i.e., with a mass matrix. These systems can represent linear DAE systems and non-proper systems.
- [TrajectoryOptimization.jl](http://roboticexplorationlab.org/TrajectoryOptimization.jl/stable/) is one of the more developed packages for open-loop **optimal control** and trajectory optimization in Julia.
- [LowLevelParticleFilters.jl](https://github.com/baggepinnen/LowLevelParticleFilters.jl) is a library for state estimation using particle filters and Kalman filters of different flavors.
- [FaultDetectionTools.jl](https://github.com/andreasvarga/FaultDetectionTools.jl) contains utilities and observers for online fault detection.
- [MatrixEquations.jl](https://github.com/andreasvarga/MatrixEquations.jl) contains solvers for many different matrix equations common in control. ControlSystems.jl makes use of this package for solving Riccati and Lyapunov equations.
- [JuMP.jl](https://jump.dev/JuMP.jl/stable/) is a modeling language for optimization, similar to YALMIP. JuMP is suitable for solving LMI/SDP problems as well as advanced linear MPC problems. 
- [SumOfSquares.jl](https://jump.dev/SumOfSquares.jl/stable/) is a package for sum-of-squares programming that builds on top of JuMP. Their documentation contains examples of Lyapunov-function search and nonlinear synthesis.
- [MonteCarloMeasurements.jl](https://baggepinnen.github.io/MonteCarloMeasurements.jl/stable/) is a library for working with parametric uncertainty. An example using ControlSystems is available [here](https://github.com/baggepinnen/MonteCarloMeasurements.jl/blob/master/examples/controlsystems.jl).
- [ModelingToolkit.jl](https://mtk.sciml.ai/stable/) is an acausal modeling tool, similar in spirit to Modelica.
- [DifferentialEquations.jl](https://diffeq.sciml.ai/stable/) is the home of the SciML ecosystem that provides solvers for scientific problems. ControlSystems.jl uses these solvers for continuous-time simulations.
- [Dojo.jl](https://github.com/dojo-sim/Dojo.jl) is a differentiable robot simulator.
- [StaticCompiler.jl](https://github.com/tshort/StaticCompiler.jl) contains tools for compiling small binaries of Julia programs.
- [JuliaPOMDP](https://github.com/JuliaPOMDP) is a Julia ecosystem for reinforcement learning. 
- [JuliaReinforcementLearning](https://github.com/JuliaReinforcementLearning) is another Julia ecosystem for reinforcement learning. 
