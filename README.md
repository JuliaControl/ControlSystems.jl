# ControlSystems.jl

[![Build Status](https://travis-ci.org/JuliaControl/ControlSystems.jl.svg?branch=master)](https://travis-ci.org/JuliaControl/ControlSystems.jl)
[![Gitter](https://badges.gitter.im/JuliaControl/ControlSystems.jl.svg)](https://gitter.im/JuliaControl/ControlSystems.jl?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge)

[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://juliacontrol.github.io/ControlSystems.jl/stable)
[![](https://img.shields.io/badge/docs-latest-blue.svg)](https://juliacontrol.github.io/ControlSystems.jl/latest)

A control systems design toolbox for Julia.

## Installation

To install, in the Julia REPL:

```julia
Pkg.add("ControlSystems")
```

## News
### 2019-05-22
New state-space type `HeteroStateSpace` that accepts matrices of heterogeneous types: [example using `StaticArrays`](https://juliacontrol.github.io/ControlSystems.jl/latest/man/creating_systems/#Creating-State-Space-Systems-1).
### 2019-01-31
System identification using [ControlSystemIdentification.jl](https://github.com/baggepinnen/ControlSystemIdentification.jl) is now available. The [readme](https://github.com/baggepinnen/ControlSystemIdentification.jl) together with a series of notebooks serve as documentation.
- [State-space identification](https://github.com/JuliaControl/ControlExamples.jl/blob/master/identification_statespace.ipynb)
- [ARX/PLR](https://github.com/JuliaControl/ControlExamples.jl/blob/master/identification_arx.ipynb)
- [Transfer-function estimation using spectral methods](https://github.com/JuliaControl/ControlExamples.jl/blob/master/identification_spectral.ipynb)
- [Impulse-response estimation](https://github.com/JuliaControl/ControlExamples.jl/blob/master/identification_impulse_response.ipynb)

### 2018-09-30
Support for Julia 0.7/1.0 added.

### 2018-09-01
- LTISystem types are now more generic and can hold matrices/vectors of arbitrary type. Examples (partly pseudo-code):
```julia
ss(1)
ss(1.)
ss(1im)
ss(ForwardDiff.Dual(1.))
ss(GPUArray([1]))
ss(SparseMatrix([1]))
```
Similar for `tf,zpk` etc.
- Continuous time systems are simulated with continuous time solvers from `OrdinaryDiffEq.jl`
- Freqresp now returns frequencies in the first dimension.
- Breaking: `lsim(sys, u::Function)` syntax has changed from `u(t,x)` to `u(x,t)` to be consistent with `OrdinaryDiffEq`
- Breaking: `feedback(P,C)` no longer returns `feedback(P*C)`. The behavior is changed to `feedback(P1, P2) = P1/(1+P1*P2)`.
- Type `Simulator` provides lower level interface to continuous time simulation.
- Example [autodiff.jl](https://github.com/JuliaControl/ControlSystems.jl/tree/master/example/autodiff.jl) provides an illustration of how the new generic types can be used for automatic differentiation of a cost function through the continuous-time solver, which allows for optimization of the cost function w.r.t. PID parameters.


## Documentation

All functions have docstrings, which can be viewed from the REPL, using for example `?tf `.

A documentation website is available at [http://juliacontrol.github.io/ControlSystems.jl/latest/](http://juliacontrol.github.io/ControlSystems.jl/latest/).

Some of the available commands are:
##### Constructing systems
ss, tf, zpk, ss2tf
##### Analysis
pole, tzero, norm, norminf, ctrb, obsv, gangoffour, margin, markovparam, damp, dampreport, zpkdata, dcgain, covar, gram, sigma, sisomargin
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
stepplot(CLs, label=["Kp = 1", "Kp = 5", "Kp = 15"])
```

![StepResponse](/example/step_response.png)

### Additional examples
See the examples folder
