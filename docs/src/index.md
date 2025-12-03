```@raw html
<p style="text-align:center">

<img src="https://avatars.githubusercontent.com/u/10605979?s=400&u=7b2efdd404c4db3b3f067f04c305d40c025a8961&v=4" alt="JuliaControl logo">

<br> 

<a class="github-button" href="https://github.com/JuliaControl/ControlSystems.jl" data-color-scheme="no-preference: light; light: light; dark: dark;" data-icon="octicon-star" data-show-count="true" aria-label="Star JuliaControl/ControlSystems.jl on GitHub">Star</a>

<script async defer src="https://buttons.github.io/buttons.js"></script>
</p> 
```

# ControlSystems.jl Manual

```@meta
CurrentModule = ControlSystems
const libpath = haskey(ENV, "CI") ? dirname(pathof(ControlSystemsBase)) : "lib/src"
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


ControlSystems.jl and the rest of the packages in the [JuliaControl](https://github.com/JuliaControl/) organization implement solutions for analysis and design of (primarily linear) control systems. If you are new to the Julia programming language, [learn more here](https://julialang.org/learning/). If you are familiar with Julia but unfamiliar with the ecosystem for control, learn more under [Ecosystem](@ref).

This documentation is structured into an introductory section labeled [Introductory guide](@ref), a section with [Examples](@ref index_examples) and a reference section sorted into topics, labeled [Functions](@ref).




## Introductory guide

```@contents
Pages = ["man/introduction.md", "man/creating_systems.md", "man/numerical.md", "man/differences.md"]
Depth = 1
```

## [Examples](@id index_examples)
```@contents
Pages = ["examples/example.md", "examples/ilc.md", "examples/smith_predictor.md"]
Depth = 2
```

## Functions

```@contents
Pages = ["lib/constructors.md",  "lib/analysis.md", "lib/synthesis.md", "lib/timefreqresponse.md", "lib/plotting.md", "lib/nonlinear.md"]
Depth = 1
```

## Ecosystem

### JuliaControl

The JuliaControl and surrounding ecosystem contains a few additional packages that may be of interest
- [RobustAndOptimalControl.jl](https://github.com/JuliaControl/RobustAndOptimalControl.jl) contains more advanced features for LQG design, robust analysis and synthesis, uncertainty modeling, named systems and an interface to [DescriptorSystems.jl](https://github.com/andreasvarga/DescriptorSystems.jl).
- [ModelPredictiveControl.jl](https://github.com/JuliaControl/ModelPredictiveControl.jl) is a package for solving both linear and nonlinear MPC problems. 
- [SymbolicControlSystems.jl](https://github.com/JuliaControl/SymbolicControlSystems.jl) contains basic C-code generation for linear systems as well as symbolic manipulation of transfer functions using SymPy.
- [ControlSystemIdentification.jl](https://github.com/baggepinnen/ControlSystemIdentification.jl) is a system-identification toolbox for identification of LTI systems using either time or frequency-domain data. This package can use data to estimate statespace models, transfer-function models and Kalman filters that can be used for control design.
- [ControlSystemsMTK.jl](https://juliacontrol.github.io/ControlSystemsMTK.jl/dev/) is an interface between ControlSystems.jl and [ModelingToolkit.jl](https://mtk.sciml.ai/stable/).
- [DiscretePIDs.jl](https://github.com/JuliaControl/DiscretePIDs.jl) contains a reference implementation in Julia of a discrete-time PID controller including set-point weighting, integrator anti-windup, derivative filtering and bumpless transfer.
- [IterativeLearningControl2.jl](https://baggepinnen.github.io/IterativeLearningControl2.jl/dev/) Implementations of several different ILC algorithms with a unified interface.

See also [the paper](https://portal.research.lu.se/en/publications/controlsystemsjl-a-control-toolbox-in-julia) describing the toolbox

> Bagge Carlson, F., Fält, M., Heimerson, A., & Troeng, O. (2021). ControlSystems.jl: A Control Toolbox in Julia. In 2021 60th IEEE Conference on Decision and Control (CDC) IEEE Press. https://doi.org/10.1109/CDC45484.2021.9683403

and the introductory Youtube video below, as well as the following [Youtube playlist](https://youtube.com/playlist?list=PLd_RaCnvGpJAk8JXI9dVtIkhoe10CadnI) with videos about using Julia for control.

```@raw html
<iframe style="height: 315px; width: 560px" src="https://www.youtube.com/embed/Fdz2Fsm1aTY" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>
```

### The wider Julia ecosystem for control
The following is a list of packages from the wider Julia ecosystem that may be of interest.

- [DescriptorSystems.jl](https://github.com/andreasvarga/DescriptorSystems.jl) contains types that represent statespace systems on descriptor form, i.e., with a mass matrix. These systems can represent linear DAE systems and non-proper systems.
- [LowLevelParticleFilters.jl](https://github.com/baggepinnen/LowLevelParticleFilters.jl) is a library for state estimation using particle filters and Kalman filters of different flavors.
- [ModelingToolkit.jl](https://mtk.sciml.ai/stable/) is an acausal modeling tool, similar in spirit to Modelica. A video showing ControlSystems and ModelingToolkit together is [available here](https://youtu.be/favQKOyyx4o). [ControlSystemsMTK.jl](https://juliacontrol.github.io/ControlSystemsMTK.jl/dev/) exists to ease the use of these two packages together.
- [DyadControlSystems.jl](https://help.juliahub.com/dyadcontrol/dev/) is a product that builds upon the JuliaControl ecosystem and ModelingToolkit, providing additional nonlinear and robust control methods.
- [LinearMPC.jl](https://darnstrom.github.io/LinearMPC.jl/stable/) A package for linear quadratic MPC that supports generation of embeddable C-code for real-time applications.
- [FaultDetectionTools.jl](https://github.com/andreasvarga/FaultDetectionTools.jl) contains utilities and observers for online fault detection.
- [ReachabilityAnalysis.jl](https://juliareach.github.io/ReachabilityAnalysis.jl/dev/generated_examples/Building/) is a package for reachability analysis. This can be used to verify stability and safety properties of linear and nonlinear systems.
- [MatrixEquations.jl](https://github.com/andreasvarga/MatrixEquations.jl) contains solvers for many different matrix equations common in control. ControlSystems.jl makes use of this package for solving Riccati and Lyapunov equations.
- [JuMP.jl](https://jump.dev/JuMP.jl/stable/) is a modeling language for optimization, similar to YALMIP. JuMP is suitable for solving LMI/SDP problems as well as advanced linear MPC problems. 
- [SumOfSquares.jl](https://jump.dev/SumOfSquares.jl/stable/) is a package for sum-of-squares programming that builds on top of JuMP. Their documentation contains examples of Lyapunov-function search and nonlinear synthesis.
- [InfiniteOpt.jl](https://github.com/infiniteopt/InfiniteOpt.jl) is a tool for solving numerical optimal-control problems built on top of JuMP.
- [MonteCarloMeasurements.jl](https://baggepinnen.github.io/MonteCarloMeasurements.jl/stable/) is a library for working with parametric uncertainty. An example using ControlSystems is available [here](https://github.com/baggepinnen/MonteCarloMeasurements.jl/blob/master/examples/controlsystems.jl).
- [DifferentialEquations.jl](https://diffeq.sciml.ai/stable/) is the home of the SciML ecosystem that provides solvers for scientific problems. ControlSystems.jl uses these solvers for continuous-time simulations.
- [StaticCompiler.jl](https://github.com/tshort/StaticCompiler.jl) contains tools for compiling small binaries of Julia programs.
- [JuliaPOMDP](https://github.com/JuliaPOMDP) is a Julia ecosystem for reinforcement learning. 
- [JuliaReinforcementLearning](https://github.com/JuliaReinforcementLearning) is another Julia ecosystem for reinforcement learning. 

## Courses using ControlSystems.jl
- [Carnegie Mellon University: Optimal-Control](https://github.com/Optimal-Control-16-745)
- [Kasetsart University: Control Engineering with Julia](https://dewdotninja.github.io/julia/control/julia_control.html)
- [Czech Technical University: Optimal and Robust Control](https://hurak.github.io/orr/). See also [github repo](https://github.com/hurak/orr/tree/main).
