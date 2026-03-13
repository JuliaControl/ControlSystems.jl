#=
This example demonstrates how to compute and plot describing functions in the Nyquist plane
using the built-in `describing_function` from ControlSystems.jl.
=#
using ControlSystems

using ControlSystemsBase: describing_function, describing_function_plot, Saturation, DeadZone, Hysteresis, nonlinearity, saturation, hysteresis


# ==============================================================================
## Examples of usage
# ==============================================================================
# 1. Sign (Relay) function: f(x) = sgn(x)
relay(x) = sign(x)
println("Relay N(2): ", describing_function(relay, 2.0))
# Theoretical check: 4/(π*A) = 4/(2π) ≈ 0.6366

# 2. Saturation function — analytical
println("Saturation N(2) (analytical): ", describing_function(Saturation(1.0), 2.0))
# Compare with numerical
println("Saturation N(2) (numerical):  ", describing_function(x -> clamp(x, -1, 1), 2.0))

# 3. Dead-zone function — analytical
println("Dead-zone N(2) (analytical): ", describing_function(DeadZone(0.1), 2.0))
# Compare with numerical
dz = DeadZone(0.1)
println("Dead-zone N(2) (numerical):  ", describing_function(x -> dz(x), 2.0))


# ==============================================================================
## Nyquist Plot with Describing Function Overlay
# ==============================================================================

using Plots

s = tf("s")
G = 10 / (s^3 + 2s^2 + s + 1)

# Create Nyquist plot with -1/N(A) overlay for saturation
describing_function_plot(G, Saturation(1.0); A_range=0.01:0.01:100)
# The intersection between the Nyquist plot and -1/N(A) indicates potential limit cycles.
# The amplitude of the limit cycle can be estimated from the corresponding A value;
# use the plotly() backend to display the amplitude on hover.


# ==============================================================================
## Simulation
# We can simulate the closed-loop system with the nonlinearity to verify the limit cycle prediction.
# ==============================================================================

nl = saturation(1.0)
sys_cl = feedback(G, nl)
plot(impulse(sys_cl, 100))


# ==============================================================================
## Hysteresis describing function
# ==============================================================================
G = 2 / ((s+1)*(s+2)) * pid(1,1,0)

H = Hysteresis(3.0, 1.5, Inf)
describing_function_plot(G, H; A_range=1.5:0.01:100)

h = hysteresis(amplitude=3.0, width=1.5, hardness=Inf)
sys_cl = feedback(G*h)
plot(impulse(sys_cl, 50))
