# In this example, we analyze the effect of ZoH sampling in continuous time and compare it to the equivalent discrete-time system
using ControlSystems, Plots
s = tf('s')
P = tf(0.1, [1, 0.1, 0.1])

Ts = 1 # Sample interval
Z = (1 - delay(Ts))/s # The transfer function of the ZoH operator
Pd = c2d(P, Ts) # The equivalent discrete-time system
Pz = P*Z # The continuous-time version of the discrete-time system

w = exp10.(-2:0.01:log10(2*0.5))

bodeplot(P, w, lab="P")
bodeplot!(Pz, w, lab="P zoh")
bodeplot!(Pd, w, lab="Pd")

vline!([0.5], l=(:black, :dash), lab="Nyquist freq.", legend=:bottomleft)

##
# The step response of Pz matches the discrete output of Pd delayed by half the sample time
resP = step(P, 12)
resPz = step(Pz, 12)
resPd = step(Pd, 12)
plot([resP, resPz, resPd], lab=["P" "Pz" "Pd"])

t_shift = resPd.t .+ Ts / 2
plot!(t_shift, resPd.y[:], lab="Pd shifted", m=:o, l=:dash)


##
# With a piecewise constant input, even if it's not a step, we get the same result
Tf = 20
ufun = (x,t)->[sin(2pi*floor(t/Ts)*Ts/5)]
resP = lsim(P, ufun, Tf)
resPz = lsim(Pz, ufun, 0:0.01:Tf)
resPd = lsim(Pd, ufun, Tf)
plot([resP, resPz, resPd], lab=["P" "Pz" "Pd"])

t_shift = resPd.t .+ Ts / 2
plot!(t_shift, resPd.y[:], lab="Pd shifted", m=:o)

##
# With a continuous input signal, the result is different,
# after the initial transient, the output of Pz matches that of Pd exactly
# try plotting with the plotly() backend and zoom in at the end
Tf = 150
ufun = (x,t)->[sin(2pi*t/5)]
resP = lsim(P, ufun, Tf)
resPz = lsim(Pz, ufun, 0:0.01:Tf)
resPd = lsim(Pd, ufun, Tf)
plot([resP, resPz, resPd], lab=["P" "Pz" "Pd"])


##
# Discrete-time systems can be converted to continuous-time systems with delays in such a way that the frequency-response of the two systems match exactly.
# To this end, the function `d2d_exact` is used. The resulting system is useful for frrequency-domain simulations only and cannot be used for time-domain simulations.
# This is useful when analyzing hybrid systems with both continuous-time and
# discrete-time components and accurate frequency-domain calculations are required over the entire frequency range up to the Nyquist frequency.
# Below, we compare the frequency response of a discrete-time system with delays to
# two continuous-tiem systems, one translated using the exact method and one using the standard `d2c` method with inverse ZOH sampling.

w = exp10.(LinRange(-2, 1, 200))
Pdc = d2c_exact(ss(Pd))
Pc = d2c(Pd)

bodeplot(Pd, w, lab="Discrete")
bodeplot!(Pdc, w, lab="Continuous with delays", l=:dash)
bodeplot!(Pc, w, lab="Continuous d2c-ZoH")
vline!([0.5 0.5], l=(:black, :dash), lab="Nyquist freq.", legend=:bottomleft)



