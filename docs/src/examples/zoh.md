# Analysis of sampled-data systems
A sampled-data system contains both continuous-tiem and discrete-time components, such as a continuous-time plant and a discrete-time controller. In this example, we will look at how to analyze such systems using the ControlSystems.jl package. To learn more about the theory of sampled-data systems, consult the reference mentioned at the end of this page.

First, we analyze the effect of ZoH sampling in continuous time and compare it to the equivalent discrete-time system

We will consider a simple second-order system ``P``, which we define in continuous time:
```@example zoh
using ControlSystems, Plots
s = tf('s')
P = tf(0.1, [1, 0.1, 0.1])
```

## Continuous to discrete

Next, we discretize this system using the standard [`c2d`](@ref) function, which uses ZoH sampling by default. We compare the frequency response of the discretized system with the frequency response of the original continuous-time system multiplied by the transfer function of the ZoH operator [^CCS]
```math
Z(s) = \dfrac{1 - e^{-T_s s}}{s}
```

```@example zoh
Ts = 1 # Sample interval
Z = (1 - delay(Ts))/s # The transfer function of the ZoH operator
Pd = c2d(P, Ts) # Discrete-time system obtained by ZoH sampling
Pz = P*Z # The continuous-time version of the discrete-time system

wd = exp10.(-2:0.01:log10(2*0.5))
bodeplot(P, wd, lab="\$P(s)\$")
bodeplot!(Pz, wd, lab="\$P(s)Z(s)\$")
bodeplot!(Pd, wd, lab="\$P_d(z)\$ (ZoH sampling)", l=:dash)
vline!([0.5 0.5], l=(:black, :dash), lab="Nyquist freq.", legend=:bottomleft)
```
The frequency response of `Pz` ``= P(s) Z(s)`` matches that of ``P_d(z)`` exactly, but these two differ from the frequency response of the original ``P(s)`` due to the ZoH operator.

The step response of `Pz` ``= P(s) Z(s)`` matches the discrete output of ``P_d(z)`` delayed by half the sample time
```@example zoh
resP = step(P, 12)
resPz = step(Pz, 12)
resPd = step(Pd, 12)
plot([resP, resPz, resPd], lab=["P" "Pz" "Pd"])

t_shift = resPd.t .+ Ts / 2
plot!(t_shift, resPd.y[:], lab="Pd shifted", m=:o, l=:dash)
```
With a piecewise constant input, even if it's not a step, we get the same result
```@example zoh
Tf = 20
ufun = (x,t)->[sin(2pi*floor(t/Ts)*Ts/5)]
resP = lsim(P, ufun, Tf)
resPz = lsim(Pz, ufun, 0:0.01:Tf)
resPd = lsim(Pd, ufun, Tf)
plot([resP, resPz, resPd], lab=["P" "Pz" "Pd"])

t_shift = resPd.t .+ Ts / 2
plot!(t_shift, resPd.y[:], lab="Pd shifted", m=:o)
```

With a _continuous_ input signal, the result is different,
after the initial transient, the output of `Pz` matches that of `Pd` exactly
(try plotting with the plotly() backend and zoom in at the end)
```@example zoh
Tf = 100
ufun = (x,t)->[sin(2pi*t/5)]
resP = lsim(P, ufun, Tf)
resPz = lsim(Pz, ufun, 0:0.01:Tf)
resPd = lsim(Pd, ufun, Tf)
plot([resP, resPz, resPd], lab=["P" "Pz" "Pd"]); lens!([Tf-10, Tf], [-0.1, 0.1], inset=(1, bbox(0.4, 0.02, 0.4, 0.3)))
```


## Discrete to continuous

Discrete-time systems can be converted to continuous-time systems formulated with delay terms in such a way that the frequency-response of the two systems match exactly, using the substitution ``z^{-1} = e^{-sT_s}``. To this end, the function [`d2c_exact`](@ref) is used. This is useful when analyzing hybrid systems with both continuous-time and discrete-time components and accurate frequency-domain calculations are required over the entire frequency range up to the Nyquist frequency.

Below, we compare the frequency response of a discrete-time system with delays to
two continuous-time systems, one translated using the exact method and one using the standard `d2c` method with inverse ZoH sampling.
We extend the frequency axis for this example to extend beyond the Nyquist frequency.
```@example zoh
wd = exp10.(LinRange(-2, 1, 200))
Pdc = d2c_exact(ss(Pd))
Pc = d2c(Pd)
bodeplot(Pd, wd, lab="\$P_d(z)\$")
bodeplot!(Pdc, wd, lab="\$P_d(s)\$ (exact translation)", l=:dash)
bodeplot!(Pc, wd, lab="\$P_d(s)\$ (inverse ZoH sampling)")
vline!([0.5 0.5], l=(:black, :dash), lab="Nyquist freq.", legend=:bottomleft)
```
We see that the translation of the discrete-time system to continuous time using the standard inverse ZoH sampling (`d2c(Pd)`) is not accurate for frequencies close to and above the Nyquist frequency. The translation using exact method (`d2c_exact(Pd)`) matches the frequency response of the discrete-time system exactly over the entire frequency axis.

## Time-domain simulation

When analyzing hybrid systems, systems with both discrete-time and continuous-time components, it is often useful to simulate the system in time-domain. To assemble the complete hybrid system, we must either
1. Convert continuous-time components to discrete time, or
2. Convert discrete-time components to continuous time.

If all inputs to the continuous-time components are piecewise constant, then the first option is the most natural. The ZoH discretization is exact for piece-wise constant inputs. This conversion can be performed using the function [`c2d`](@ref) with the (default) ZoH sampling method. If some of the inputs to the continuous-time components are continuously varying, then the second option may be used for maximum accuracy. This is particularly useful if the continuous-time inputs contain frequencies higher than the Nyquist frequency, or when the continuous-time plant contains significant dynamics (resonances etc.) at frequencies higher than the Nyquist frequency of the discrete-time part. This conversion can be performed using the function [`d2c_exact`](@ref) which has two modes of operation. The default method produces a causal system which can be simulated in the time domain, while the second method produces an acausal system of lower complexity, which is amenable to frequency-domain computations only. Since we are going to simulate in the time domain here, we will use the causal method (default).

In this example, we will model a continuous-time system ``P(s)`` with resonances appearing at large frequencies, while the discrete-time control system is operating a significantly lower frequencies.
```@example zoh
A = [0.0    0.0    0.0   1.0     0.0     0.0
   0.0    0.0    0.0   0.0     1.0     0.0
   0.0    0.0    0.0   0.0     0.0     1.0
 -10.0    0.0   10.0  -0.001   0.0     0.001
  -0.0  -10.0   10.0  -0.0    -0.001   0.001
  10.0   10.0  -20.0   0.001   0.001  -0.002]
B = [0, 0, 0, 0.1, 0, 0]
C = [0.0  0.0  0.0  0.0  0.0  1.0]
P = minreal(ss(A,B,C,0))
w = exp10.(LinRange(-2, 2, 1000))
bodeplot(P, w, lab="\$P(s)\$", plotphase=false)
```
The discrete-time controller ``C(z)`` is a basic PI controller

```@example zoh
Ts = 1
C = pid(0.01,10; Ts, Tf = 1/100, state_space=true)
```
If we choose option 1. in this case, we get a discretized plant that has a very poor frequency-response fit to the original continuous-time plant, making frequency-domain analysis difficult
```@example zoh
Pd = c2d(P, Ts)
Ld = Pd*C
wd = exp10.(-3:0.001:log10(1000*0.05))
figb = bodeplot(P, wd, lab="\$P(s)\$")
bodeplot!(Pd, wd, lab="\$P(z)\$ ZoH")
bodeplot!(Ld, wd, lab="\$P(z)C(z)\$", legend=:bottomleft)
```

If we instead make use of the second method above, exactly translating the discrete-time controller to continuous time, we can more easily determine that the closed-loop system will be stable by inspecting the Bode plot. Here, we show the Bode plot of ``P`` and that of the loop-transfer function, including the ZoH operator, ``P(s)Z(s)C(s)``.
```@example zoh
Cc = d2c_exact(C)
Lc = P*Z*Cc

bodeplot(P, wd, lab="\$P(s)\$")
bodeplot!(Lc, wd, lab="\$P(s)C(s)\$", legend=:bottomleft, c=3)
```

If we form the closed-loop system from an input disturbance to the output
```math
PS = \dfrac{P(s)}{1 + P(s)C(s)}
```
we can simulate the response to a continuous-time input disturbance that contains frequencies higher than the Nyquist frequency of the discrete-time system, we do this below. We also try doing this with the discretized plant, which yields a very poor result
```@example zoh
PS = feedback(P, Z*Cc)  # Use the continuous-time plant and continuous translation of the controller + ZoH operator
PSd = feedback(Pd, C) # Use the discretized plant and discrete controller
ω = 5.53 # Frequency of input disturbance (rad/s) (Nyquist frequency is π rad/s)
disturbance(x, t) = sin(ω*t) # Continuous input disturbance
plot(lsim(PS, disturbance, 0:0.22:3500), lab="Continuous disturbance response")
plot!(lsim(PSd, disturbance, 3500), lab="Discrete disturbance response")
hline!([abs(freqresp(PS, ω)[])], l=(:black, :dash), lab="Predicted freq. response amplitude", legend=:bottomleft, fmt=:png)
```
The continuous-time analysis eventually settles at a periodic output that matches the amplitude predicted by the continuous-time frequency response. However, the discrete-time simulation yields, as expected, a very poor result. Noticeable in the simulation is the appearance of a beat frequency, which is expected due to the modulation introduced by sampling.[^CCS]

## Caveats
- The exact output of the system that was translated from discrete to continuous time is not going to be accurate, so transient properties of the hybrid system cannot be accurately assessed using this kind of analysis. 
- Interpretation of frequency-responses for sampled-data systems must be done with care. The modulation introduced by sampling implies the creating of additional frequencies in the output. For an input with frequency ``\omega``, the output will contain all frequencies ``\omega ± \omega_s k`` where ``\omega_s`` is the sampling frequency and ``k`` is an integer.[^CCS]

## References
Learn more about sampled-data systems and zero-order hold sampling in chapter 7 of the reference below.
[^CCS]: Åström, K. J., & Wittenmark, B. (1997). Computer-Controlled Systems: Theory and Design.

