# ControlSystems.jl

[![Build Status](https://travis-ci.org/JuliaControl/ControlSystems.jl.svg?branch=master)](https://travis-ci.org/JuliaControl/ControlSystems.jl)
[![Coverage Status](https://coveralls.io/repos/github/JuliaControl/ControlSystems.jl/badge.svg?branch=master)](https://coveralls.io/github/JuliaControl/ControlSystems.jl?branch=master)

A control systems design toolbox for Julia.

## Installation

To install, in the Julia REPL:

```julia
Pkg.add("ControlSystems")
```

Note that the latest version of this package requires Julia 0.5. Users of Julia 0.4 should use [v0.1.4](https://github.com/JuliaControl/ControlSystems.jl/tree/v0.1.4) instead.

## Documentation

All functions have docstrings, which can be viewed from the REPL, using for example `?tf `.

A documentation website under developement is available at [http://juliacontrol.github.io/ControlSystems.jl/latest/](http://juliacontrol.github.io/ControlSystems.jl/latest/).

Some of the available commands are:
##### Constructing systems
ss, tf, zpk, ss2tf
##### Analysis
pole, tzero, norm, norminf, ctrb, obsv, gangoffour, margin, markovparam, damp, dampreport, zpkdata, dcgain, covar, gram, sigma, sisomargin
##### Synthesis
care, dare, dlyap, lqr, dlqr, place, pid, leadlink, laglink, leadlinkat, rstd, rstc, dab
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
stepplot(CLs);
# Add legend. This is specific to the PyPlot back-end, see documentation for more info.
PyPlot.legend(["Kp = 1", "Kp = 5", "Kp = 15"]);
```

![StepResponse](/example/step_response.png)
