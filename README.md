# ControlSystems.jl

[![Build Status](https://github.com/JuliaControl/ControlSystems.jl/workflows/CI/badge.svg)](https://github.com/JuliaControl/ControlSystems.jl/actions?query=workflow%3ACI)
[![Documentation Status](https://github.com/JuliaControl/ControlSystems.jl/workflows/Docs/badge.svg)](https://github.com/JuliaControl/ControlSystems.jl/actions?query=workflow%3ADocs)

[![PkgEval](https://juliaci.github.io/NanosoldierReports/pkgeval_badges/C/ControlSystems.svg)](https://juliaci.github.io/NanosoldierReports/pkgeval_badges/report.html)
[![codecov](https://codecov.io/gh/JuliaControl/ControlSystems.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/JuliaControl/ControlSystems.jl)

[![](https://img.shields.io/badge/docs-latest-blue.svg)](https://juliacontrol.github.io/ControlSystems.jl/latest)
[![Star on GitHub](https://img.shields.io/github/stars/JuliaControl/ControlSystems.jl.svg?style=social)](https://github.com/JuliaControl/ControlSystems.jl/stargazers)

A control systems design toolbox for Julia.

## Installation

To install, in the Julia REPL:

```julia
using Pkg; Pkg.add("ControlSystems")
```

## News

### 2022-09 v1.5
ControlSystems have now been divided into two packages:
- ControlSystemsBase.jl lives in this same repository in the folder `lib/`, and contains most of the original functionality of ControlSystems, except continuous-time simulation and root locus.
- ControlSystems.jl internally uses and reexports ControlSystemsBase.jl and adds continuous-time simulation.

This is **not** a breaking change, users of ControlSystems.jl will have the same functionality now as before, but users who do not need continuous-time simulation (including simulation of delay systems and nonlinear systems) can use the considerably more lightweight ControlSystemsBase.jl package instead. 

OrdinaryDiffEq.jl and DelayDiffEq.jl contributed the vast majority of both pre-compilation time and loading time of ControlSystems.jl, and workflows that do not require these packages were thus burdened by this overhead unnecessarily. If you do not need this, install and use ControlSystemsBase rather than ControlSystems (do not even install ControlSystems to avoid the pre-compilation time of OrdinaryDiffEq).

We still encourage advanced users to build a system image following [the instructions in the documentation](https://juliacontrol.github.io/ControlSystems.jl/latest/man/differences/#Precompilation-for-faster-load-times), but in the absence of such a system image, we now have the following timings:
```julia
julia> @time using ControlSystemsBase
  1.70 seconds
```
```julia
julia> @time using ControlSystems
  10.90 seconds
```


### 2022-09 v1.4
New feature: `loopshapingPID`
Release notes: https://github.com/JuliaControl/ControlSystems.jl/releases

### 2022-07 v1.2
Better support for static systems (using StaticArrays)

### 2022-07 v1.0
- *Breaking*: Frequency-responses have changed data layout to `ny×nu×nω` from the previous `nω×ny×nu`. This is for performance reasons and to be consistent with time responses. This affects downstream functions `bode` and `nyquist` as well.
- *Breaking*: `baltrunc` and `balreal` now return the diagonal of the Gramian as the second argument rather than the full matrix.
- *Breaking*: The `pid` constructor no longer takes parameters as keyword arguments. `pid` has also gotten some new features, the new signature is `pid(P, I, D=0; form = :standard, Ts=nothing, Tf=nothing, state_space=false)`. This change affects downstream functions like `placePI, loopshapingPI, pidplots`.
- *Breaking*: The semantics of broadcasted multiplication between two systems was previously inconsistent between `StateSpace` and `TransferFunction`. The new behavior is documented under [Multiplying systems](https://juliacontrol.github.io/ControlSystems.jl/latest/man/creating_systems/#Multiplying-systems) in the documentation.

## Documentation

All functions have docstrings, which can be viewed from the REPL, using for example `?tf `.

A documentation website is available at [http://juliacontrol.github.io/ControlSystems.jl/latest/](http://juliacontrol.github.io/ControlSystems.jl/latest/) and an [introductory video is available here](https://www.youtube.com/watch?v=Fdz2Fsm1aTY&ab_channel=jolin%E2%80%A4io).

Some of the available commands are:
##### Constructing systems
`ss, tf, zpk, delay`
##### Analysis
`poles, tzeros, norm, hinfnorm, linfnorm, ctrb, obsv, gangoffour, margin, markovparam, damp, dampreport, zpkdata, dcgain, covar, gram, sigma, sisomargin`
##### Synthesis
`are, lyap, lqr, place, leadlink, laglink, leadlinkat, rstd, rstc, dab, balreal, baltrunc`
###### PID design
`pid, stabregionPID, loopshapingPI, pidplots, placePI`
##### Time and Frequency response
`step, impulse, lsim, freqresp, evalfr, bode, nyquist`
##### Plotting
`bodeplot, nyquistplot, sigmaplot, plot(lsim(...)), plot(step(...)), plot(impulse(...)), marginplot, gangoffourplot, pzmap, nicholsplot, pidplots, rlocus, leadlinkcurve`
##### Other
`minreal, sminreal, c2d`
## Usage

This toolbox works similar to that of other major computer-aided control
systems design (CACSD) toolboxes. Systems can be created in either a [transfer
function](http://en.wikipedia.org/wiki/Transfer_function) or a [state
space](http://en.wikipedia.org/wiki/State-space_representation) representation.
These systems can then be combined into larger architectures, simulated in both
time and frequency domain, and analyzed for stability/performance properties.

### Example

Here we create a simple position controller for an electric motor with an
inertial load.

```julia
using ControlSystems

# Motor parameters
J = 2.0
b = 0.04
K = 1.0
R = 0.08
L = 1e-4

# Create the model transfer function
s = tf("s")
P = K/(s*((J*s + b)*(L*s + R) + K^2))
# This generates the system
# TransferFunction:
#                1.0
# ---------------------------------
# 0.0002s^3 + 0.160004s^2 + 1.0032s
#
#Continuous-time transfer function model

# Create an array of closed loop systems for different values of Kp
CLs = TransferFunction[kp*P/(1 + kp*P) for kp = [1, 5, 15]];

# Plot the step response of the controllers
# Any keyword arguments supported in Plots.jl can be supplied
using Plots
plot(step.(CLs, 5), label=["Kp = 1" "Kp = 5" "Kp = 15"])
```

![StepResponse](/example/step_response.png)

### Additional examples
See the examples folder and [ControlExamples.jl](https://github.com/JuliaControl/ControlExamples.jl/) and several examples in the [documentation](http://juliacontrol.github.io/ControlSystems.jl/latest/).
