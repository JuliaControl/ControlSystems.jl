# Properties of delay systems
Delay systems can sometimes have non-intuitive properties, in particular when the delays appear inside of the system, i.e., not directly on the inputs or outputs. 

The Nyquist plot of delay systems usually spirals towards the origin for delay systems. This is due to the phase loss at high frequencies due to the delay:
```@example DELAY
using ControlSystemsBase, Plots
w = exp10.(LinRange(-2, 2, 2000))
P = tf(1, [1, 1]) * delay(2) # Plant with delay on the input
nyquistplot(P, w)
```

When forming a feedback interconnection, making the delay appear in the closed loop, we may get gain ripple:
```@example DELAY
bodeplot(feedback(P), w)
```

If the system with delay has a direct feedthrough term, step responses may show repeated steps at integer multiples of the delay:
```@example DELAY
using ControlSystems # Load full control systems to get simulation functionality
P = tf([1, 1], [1, 0])*delay(1)
plot(step(feedback(P, 0.5), 0:0.001:20))
```
Indeed, if the system has a non-zero feedthrough, the output will contain a delayed step attenuated by the feedthrough term, in this case
```@example DELAY
ss(feedback(tf([1, 1], [1, 0]))).D[]
```
the steps will thus in this case decay exponentially with decay rate 0.5. 

For a more advanced example using time delays, see the [Smith predictor](@ref) tutorial.

## Simulation of time-delay systems
Time-delay systems are numerically challenging to simulate, if you run into problems, please open an issue with a reproducing example. The [`lsim`](@ref), [`step`](@ref) and [`impulse`](@ref) functions accept keyword arguments that are passed along to the ODE integrator, this can be used to both select integration method and to tweak the integrator options. The documentation for solving delay-differential equations is available [here](https://diffeq.sciml.ai/latest/types/dde_types/) and [here](https://diffeq.sciml.ai/latest/tutorials/dde_example/).

## Estimation of delay
See the companion tutorial in ControlSystemIdentification.jl on [Delay estimation](https://baggepinnen.github.io/ControlSystemIdentification.jl/dev/examples/delayest/). This tutorial covers the both the detection of the presence of a delay, and estimation of models for systems with delays.

## Approximation and discretization of delays
Delay systems may be approximated as rational functions by means of [Padé approximation](https://en.wikipedia.org/wiki/Pad%C3%A9_approximant) using the function [`pade`](@ref). Pure continuous-time delays can also be discretized using the function [`thiran`](@ref). Continuous-time models with internal delays can be discretized using [`c2d`](@ref), provided that the delay is an integer multiple of the sampling time (fractional delays are not yet supported by [`c2d`](@ref)).