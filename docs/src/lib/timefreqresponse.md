# Time and Frequency response analysis
```@index
Pages = ["timefreqresponse.md"]
```



## Frequency response

Frequency responses are calculated using [`freqresp`](@ref), [`bode`](@ref), [`sigma`](@ref) and [`nyquist`](@ref). Frequency-response plots are obtained using [`bodeplot`](@ref), [`sigmaplot`](@ref), [`nyquistplot`](@ref), [`marginplot`](@ref) and [`nicholsplot`](@ref).

Any [`TransferFunction`](@ref) can be evaluated at a point using
`F(s)`, `F(omega, true)`, `F(z, false)`

- `F(s)` evaluates the continuous-time transfer function `F` at `s`.
- `F(omega,true)` evaluates the discrete-time transfer function `F` at `exp(i*Ts*omega)`
- `F(z,false)` evaluates the discrete-time transfer function `F` at `z`

A video demonstrating frequency-response analysis in ControlSystems.jl is available below.

```@raw html
<iframe style="height: 315px; width: 560px" src="https://www.youtube.com/embed/_ZvvRLpCLG0" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>
```

## Time response (simulation)

Simulation with arbitrary inputs is primarily handled by the function [`lsim`](@ref), with [`step`](@ref) and [`impulse`](@ref) serving as convenience functions to simulate responses to particular inputs.

The function [`lsim`](@ref) can take an input vector `u` containing a sampled input trajectory, or an input function taking the state and time as arguments, `u(x,t)`. This function can be used to easily simulate, e.g., ramp responses or saturated state-feedback control etc. See the docstring of [`lsim`](@ref) for more details.

For more extensive nonlinear simulation capabilities, see the notes on ModelingToolkit and DifferentialEquations under [The wider Julia ecosystem for control](@ref).

#### Example:
The following simulates a step response of a second-order system and plots the result.
```@example TIMERESP
using ControlSystemsBase, Plots
G = tf(1, [1, 1, 1])
res = step(G, 20) # Simulate 20 seconds step response
plot(res)
```
Using the function [`stepinfo`](@ref), we can compute characteristics of a step response:
```@example TIMERESP
si = stepinfo(res)
```

We can also plot the [`StepInfo`](@ref) object
```@example TIMERESP
plot(si)
```

A video demonstrating time-domain simulation in ControlSystems.jl is available below.
```@raw html
<iframe style="height: 315px; width: 560px" src="https://www.youtube.com/embed/ii86sg_8xGw" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>
```

## Docstrings

```@autodocs
Modules = [ControlSystems, ControlSystemsBase]
Pages   = [libpath*"/timeresp.jl", libpath*"/result_types.jl", libpath*"/freqresp.jl", "simulators.jl"]
Order   = [:function, :type]
Private = false
```
