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

#### Example step response:
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

#### Example `lsim`:

The function [`lsim`](@ref) can take the control input as either
1. An array of equidistantly sampled values, in this case the argument `u` is expected to have the shape `nu × n_time`
2. A function of the state and time `u(x,t)`. This form allows simulation of state feedback, a step response at time ``t_0``: `u(x, t) = amplitude * (t > t0)`, or a ramp response: `u(x, t) = t` etc.

The example below simulates state feedback with a step disturbance at ``t=4`` by providing the function `u(x,t) = -L*x .+ (t > 4)` to `lsim`:
```@example TIMERESP
using ControlSystems
using LinearAlgebra: I
using Plots

A = [0 1; 0 0]
B = [0;1]
C = [1 0]
sys = ss(A,B,C,0)
Q = I
R = I
L = lqr(sys,Q,R)

u(x,t) = -L*x .+ (t > 4) # State feedback + step disturbance
t  = 0:0.1:12
x0 = [1,0]
y, t, x, uout = lsim(sys,u,t,x0=x0)
plot(t,x', lab=["Position" "Velocity"], xlabel="Time [s]"); vline!([4], lab="Step disturbance", l=(:black, :dash, 0.5))
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
