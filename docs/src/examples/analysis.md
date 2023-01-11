# Analysis of linear control systems
From classical control, we get robustness measures such as gain and phase margins. These provide a quick and intuitive way to assess robustness of single-input, single-output systems, but also have a number of downsides, such as optimism in the presence of simultaneous gain and phase variations as well as limited applicability for MIMO systems.

Gain and phase margins can be computed using the functions [`margin`](@ref) and [`marginplot`](@ref)

## Example: Gain and phase margins
```@example
using ControlSystemsBase, Plots
P = tf(1, [1, 0.2, 1])
C = pid(0.2, 1)
loopgain = P*C
marginplot(loopgain)
```
This plot tells us that there is one gain margin of 1.27, i.e., the gain can increase by a factor of 1.27 before the system goes unstable. It also tells us that there are three different phase margins, the smallest of which is about 9°. We usually aim for a gain margin of >1.5 and a phase margin above 30-45° for a robust system. The vertical lines in the plot indicate the frequencies at which the margins have been computed.

## Sensitivity analysis
More generally applicable measures of robustness include analysis of sensitivity functions, notably the peaks of the sensitivity function
```math
S(s) = (I + P(s)C(s))^{-1}
```
and the complementary sensitivity function
```math
T(s) = I - S(s) = (I + P(s)C(s))^{-1}P(s)C(s)
```

### Examples
We can plot all four sensitivity functions referred to as the "gang of four" using [`gangoffourplot`](@ref).
```@example SENS
using ControlSystemsBase, Plots
P = tf(1, [1, 0.2, 1])
C = pid(0.2, 1)
gangoffourplot(P, C)
```

The peak value of the sensitivity function, ``M_S``, can be computed using [`hinfnorm`](@ref)
```@example SENS
S = sensitivity(P, C)
Ms, ωMs = hinfnorm(S)
```

And we can plot a circle in the Nyquist plot corresponding to the inverse distance between the loop-transfer function and the critical point:
```@example SENS
w = exp10.(-1:0.001:2)
nyquistplot(P*C, w, Ms_circles=[Ms], xlims=(-1.2, 0.5), ylims=(-2, 0.3))
```

``M_S`` is always ``≥ 1``, but we typically want to keep it below 1.3-2 for robustness reasons. For SISO systems, ``M_S`` is linked to the classical gain and phase margins through the following inequalities:
```math
\begin{aligned}
\phi_m &≥ 2 \sin^{-1}\left(\dfrac{1}{2M_S}\right) \text{rad}\\
g_m &≥ \dfrac{M_S}{M_S-1}
\end{aligned}
```

We can also obtain individual sensitivity function using the low-level function [`feedback`](@ref) directly, or using one of the higher-level functions
- [`sensitivity`](@ref)
- [`comp_sensitivity`](@ref)
- [`G_PS`](@ref)
- [`G_CS`](@ref)
- [`gangoffour`](@ref)
- [`extended_gangoffour`](@ref)
- [`feedback_control`](@ref)


## Further reading
A modern robustness measure is the [`diskmargin`](https://juliacontrol.github.io/RobustAndOptimalControl.jl/dev/#Diskmargin-example), that analyses the robustness of a SISO or MIMO system to simultaneous gain and phase variations.

In the presence of structured uncertainty, such as parameter uncertainty or other explicitly modeled uncertainty, the structured singular value (often referred to as $\mu$), provides a way to analyze robustness with respect to the modeled uncertainty. See the [RobustAndOptimalControl.jl](https://juliacontrol.github.io/RobustAndOptimalControl.jl/dev/) package for more details.