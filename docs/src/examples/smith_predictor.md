# Smith predictor
This example designs a controller for a plant with a time delay using the internal-model principle, which in this case implies the use of a Smith predictor. The plant is given by
$$ \dfrac{1}{s + 1}e^{-s\tau} = P_0 e^{-s\tau}$$

and the control architecture looks like this
```
                ┌──────┐              ┌─────────────┐
r               │      │          u   │             │
───+──+────────►│  C0  ├───────────┬─►│ P0*exp(-st) ├─┐y
   ▲  ▲         │      │           │  │             │ │
  -│  │-        └──────┘           │  └─────────────┘ │
   │  │                            │                  │
   │  │ ┌──────────┐    ┌──────┐   │                  │
   │  │ │          │    │      │   │                  │
   │  └─┤1-exp(-st)│◄───┤  P0  │◄──┘                  │
   │    │          │    │      │                      │
   │    └──────────┘    └──────┘                      │
   │                                                  │
   └──────────────────────────────────────────────────┘
```
The benefit of this approach is that the controller $C_0$ can be designed for the nominal plant $P_0$ without time delay, and still behave well in the presence of the delay. We also see why we refer to such a controller as using an "internal model", due to the presence of a model of $P_0$ in the inner feedback path.

We now set up the nominal system and PI controller

```@example smith
using ControlSystemsBase, Plots
P0 = ss(-1, 1, 1, 0) # Nominal system
```

We design a PI controller for nominal system using [`placePI`](@ref). To verify the pole placement, use, e.g., `dampreport(feedback(P0, C0))`
```@example smith
ω0 = 2
ζ  = 0.7
C0, _ = placePI(P0, ω0, ζ)
```
We then setup delayed plant + Smith predictor-based controller
```@example smith
τ = 8
P = delay(τ) * P0
C = feedback(C0, (1.0 - delay(τ))*P0) # form the inner feedback connection in the diagram above
```
We now plot the closed loop responses. The transfer function from $r$ to $y$ is given by $PC_r/(1+PC_r)$ = `feedback(P*C,1)`, and from a load disturbance entering at $u$ the transfer function is $P/(1+PC_r)$ = `feedback(P, C)`
```@example smith
using ControlSystems # Load full ControlSystems for delay-system simulation
G = [feedback(P*C, 1) feedback(P, C)] # Reference step at t = 0 and load disturbance step at t = 15
fig_timeresp = plot(lsim(G, (_,t) -> [1; t >= 15], 0:0.1:40),  title="τ = $τ")
```
Plot the frequency response of the predictor part and compare to a negative delay, which would be an ideal controller that can (typically) not be realized in practice (a negative delay implies foresight). 
```@example smith
C_pred = feedback(1, C0*(ss(1.0) - delay(τ))*P0)
fig_bode = bodeplot([C_pred, delay(-τ)], exp10.(-1:0.002:0.4), ls=[:solid :solid :dash :dash], title="", lab=["Smith predictor" "" "Ideal predictor" ""])
plot!(yticks=[0.1, 1, 10], sp=1)
plot!(yticks=0:180:1080, sp=2)
```
Check the Nyquist plot. Note that the Nyquist curve encircles -1 for τ > 2.99
```@example smith
fig_nyquist = nyquistplot(C * P, exp10.(-1:1e-4:2), title="τ = $τ")
```

A video tutorial on delay systems is available here:
```@raw html
<iframe style="height: 315px; width: 560px" src="https://www.youtube.com/embed/ER8_oHU2vZs" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>
```

## Additional design methods for delay systems
Many standard control-design methods fail for delay systems, or any system not represented as a rational function. In addition to using the Smith predictor outlined above, there are however several common tricks that can be applied to make use of these methods.
- Approximate the delay using a [`pade`](@ref) approximation, this will result in a standard rational model. The drawbacks include zeros in the right half plane and a failure to capture the extreme phase loss of the delay for high frequencies.
- Discretize the system with a sample time that fits an integer multiple in the delay time. A delay can be represented exactly in discrete time, but if the sample time is chosen small in relation to the delay time, a large number of extra states will be introduced.
- Neglect the delay and design the controller with large phase and delay margins. This is perhaps not a terribly sophisticated method, but nevertheless useful in practice.
- Neglect the delay, but model it as uncertainty. See [Modeling uncertain time delays](https://juliacontrol.github.io/RobustAndOptimalControl.jl/dev/uncertainty/#Uncertain-time-delays) in the RobustAndOptimalControl.jl extension package. This can help you get a feeling for the margin with which you must design your controller when you have neglected to model the delay.
- Frequency-domain methods such as manual loop shaping, and some forms of optimization-based tuning, handle time delays natively. 

Whatever method is used to design in the presence of delays, the robustness and performance of the design should preferably be verified using a model of the plant where the delay is included, uncertain or not.