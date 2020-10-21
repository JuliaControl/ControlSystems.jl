# ControlSystems.jl

[![Build Status](https://travis-ci.com/JuliaControl/ControlSystems.jl.svg?branch=master)](https://travis-ci.com/JuliaControl/ControlSystems.jl)
[![PkgEval](https://juliaci.github.io/NanosoldierReports/pkgeval_badges/C/ControlSystems.svg)](https://juliaci.github.io/NanosoldierReports/pkgeval_badges/report.html)
[![codecov](https://codecov.io/gh/JuliaControl/ControlSystems.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/JuliaControl/ControlSystems.jl)

[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://juliacontrol.github.io/ControlSystems.jl/stable)
[![](https://img.shields.io/badge/docs-latest-blue.svg)](https://juliacontrol.github.io/ControlSystems.jl/latest)

A control systems design toolbox for Julia.

## Installation

To install, in the Julia REPL:

```julia
using Pkg; Pkg.add("ControlSystems")
```

## News

### 2020-10
- `lsimplot, stepplot, impulseplot` now have the same signatures as the corresponding non-plotting function.
- New function `d2c` for conversion from discrete to continuous.

### 2020-09-24
Release v0.7 introduces a new `TimeEvolution` type to handle `Discrete/Continuous` systems. See the [release notes](https://github.com/JuliaControl/ControlSystems.jl/releases/tag/v0.7.0).

### 2019-11-03
- Poles and zeros are "not sorted" as in Julia versions < 1.2, even on newer versions of Julia. This should imply that complex conjugates are kept together.

### 2019-05-28
#### Delay systems
- We now support systems with time delays. Example:
```julia
sys = tf(1, [1,1])*delay(1)
stepplot(sys, 5) # Compilation time might be long for first simulation
nyquistplot(sys)
```
#### New examples
- [Delayed systems (frequency domain)](https://github.com/JuliaControl/ControlSystems.jl/blob/master/example/delayed_lti_system.jl)
- [Delayed systems (time domain)](https://github.com/JuliaControl/ControlSystems.jl/blob/master/example/delayed_lti_timeresp.jl)
- [Systems with uncertainty](https://github.com/baggepinnen/MonteCarloMeasurements.jl/blob/master/examples/controlsystems.jl)
- [Robust PID optimization](https://github.com/baggepinnen/MonteCarloMeasurements.jl/blob/master/examples/robust_controller_opt.jl)
### 2019-05-22
New state-space type `HeteroStateSpace` that accepts matrices of heterogeneous types: [example using `StaticArrays`](https://juliacontrol.github.io/ControlSystems.jl/latest/man/creating_systems/#Creating-State-Space-Systems-1).
### 2019-01-31
System identification using [ControlSystemIdentification.jl](https://github.com/baggepinnen/ControlSystemIdentification.jl) is now available. The [readme](https://github.com/baggepinnen/ControlSystemIdentification.jl) together with a series of notebooks serve as documentation.
- [State-space identification](https://github.com/JuliaControl/ControlExamples.jl/blob/master/identification_statespace.ipynb)
- [ARX/PLR](https://github.com/JuliaControl/ControlExamples.jl/blob/master/identification_arx.ipynb)
- [Transfer-function estimation using spectral methods](https://github.com/JuliaControl/ControlExamples.jl/blob/master/identification_spectral.ipynb)
- [Impulse-response estimation](https://github.com/JuliaControl/ControlExamples.jl/blob/master/identification_impulse_response.ipynb)


## Documentation

All functions have docstrings, which can be viewed from the REPL, using for example `?tf `.

A documentation website is available at [http://juliacontrol.github.io/ControlSystems.jl/latest/](http://juliacontrol.github.io/ControlSystems.jl/latest/).

Some of the available commands are:
##### Constructing systems
ss, tf, zpk
##### Analysis
pole, tzero, norm, hinfnorm, linfnorm, ctrb, obsv, gangoffour, margin, markovparam, damp, dampreport, zpkdata, dcgain, covar, gram, sigma, sisomargin
##### Synthesis
care, dare, dlyap, lqr, dlqr, place, leadlink, laglink, leadlinkat, rstd, rstc, dab, balreal, baltrunc
###### PID design
pid, stabregionPID, loopshapingPI, pidplots
##### Time and Frequency response
step, impulse, lsim, freqresp, evalfr, bode, nyquist
##### Plotting
lsimplot, stepplot, impulseplot, bodeplot, nyquistplot, sigmaplot, marginplot, gangoffourplot, pidplots, pzmap, nicholsplot, pidplots, rlocus, leadlinkcurve
##### Other
minreal, sminreal, c2d
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
stepplot(CLs, label=["Kp = 1" "Kp = 5" "Kp = 15"])
```

![StepResponse](/example/step_response.png)

### Additional examples
See the examples folder
