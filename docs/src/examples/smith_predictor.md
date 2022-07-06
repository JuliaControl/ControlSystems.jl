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
using ControlSystems, Plots
P0 = ss(-1, 1, 1, 0) # Nominal system
```

We design a PI controller for nominal system using [`placePI`](@ref). To verify the pole placement, use, e.g., `dampreport(feedback(P0, C0))`
```@example smith
ω0 = 2
ζ  = 0.7
_, C0 = placePI(P0, ω0, ζ)
```
We then setup delayed plant + Smith predictor-based controller
```@example smith
τ = 8
P = delay(τ) * P0
C = feedback(C0, (1.0 - delay(τ))*P0) # form the inner feedback connection in the diagram above
```
We now plot the closed loop responses. The transfer function from $r$ to $y$ is given by $PC_r/(1+PC_r)$ = `feedback(P*C,1)`, and from a load disturbance entering at $u$ the transfer function is $P/(1+PC_r)$ = `feedback(P, C)`
```@example smith
G = [feedback(P*C, 1) feedback(P, C)] # Reference step at t = 0 and load disturbance step at t = 15
fig_timeresp = plot(lsim(G, t -> [1; t >= 15], 0:0.1:40),  title="τ = $τ")
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
