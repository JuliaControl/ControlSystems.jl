# Control.jl

[![Build Status](https://travis-ci.org/JuliaControl/Control.jl.svg)](https://travis-ci.org/JuliaControl/Control.jl)

A control systems design toolbox for Julia.

## Installation

This package is yet to be released, and as such is not in `METADATA` (it should
be released soon though!). To install, in the Julia REPL:

```julia
julia> Pkg.clone("https://github.com/JuliaControl/Control.jl.git")
```

Note that this package requires Julia 0.4.

## Documentation

All exported functions have docstrings, which can be viewed from the REPL. A
documentation website is on the list of things to do...

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
julia> using Control

# Motor parameters
julia> J = 2.0

julia> b = 0.04

julia> K = 1.0

julia> R = 0.08

julia> L = 1e-4

# Create the model transfer function
julia> s = tf("s")

julia> P = K/(s*((J*s + b)*(L*s + R) + K^2))
TransferFunction:
               1.0
---------------------------------
0.0002s^3 + 0.160004s^2 + 1.0032s

Continuous-time transfer function model

# Create an array of closed loop systems for different values of Kp
julia> CLs = LTISystem[kp*P/(1 + kp*P) for kp = [1, 5, 15]];

# Plot the step response of the controllers
julia> stepplot(CLs);

julia> legend(["Kp = 1", "Kp = 5", "Kp = 15"]);
```

![StepResponse](/example/step_response.png)
