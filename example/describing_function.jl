#=
This example demonstrates how to compute and plot describing functions in the Nyquist plane
=#
using ControlSystems
using QuadGK

"""
    compute_describing_function(f, A)

Numerically computes the describing function N(A) for a nonlinearity f(x)
at input amplitude A. Returns a Complex number.
"""
function compute_describing_function(f, A)
    # Fundamental sine component coefficient (b1)
    # b1 = (1/π) * ∫ [f(A*sin(θ)) * sin(θ)] dθ from 0 to 2π
    b1_integral, _ = quadgk(θ -> f(A * sin(θ)) * sin(θ), 0, 2π, atol=1e-6, rtol=1e-6)
    b1 = b1_integral / π

    # Fundamental cosine component coefficient (a1)
    # a1 = (1/π) * ∫ [f(A*sin(θ)) * cos(θ)] dθ from 0 to 2π
    a1_integral, _ = quadgk(θ -> f(A * sin(θ)) * cos(θ), 0, 2π, atol=1e-6, rtol=1e-6)
    a1 = a1_integral / π

    # N(A) = (b1 + j*a1) / A
    return Complex(b1, a1) / A
end

"""
    plot_negative_inverse_df(f; A_range=range(0.1, 10, length=100), label="Nonlinearity")
Computes and plots -1/N(A) for a given function f over a range of amplitudes.
"""
function plot_negative_inverse_df(f; A_range=range(0.1, 10, length=100), label="Nonlinearity", fig=nothing)
    # Compute the negative inverse values
    vals = similar(A_range, ComplexF64)
    Threads.@threads for i = eachindex(A_range)
        A = A_range[i]
        vals[i] = -1.0 / compute_describing_function(f, A)
    end
    
    # Extract real and imaginary parts for plotting
    re_vals = real.(vals)
    im_vals = imag.(vals)
    
    # Create plot or add to existing one
    if fig === nothing
        fig = plot(title="Negative Inverse Describing Function Analysis", 
                   xlabel="Real", ylabel="Imaginary", grid=true, zeroline=true)
    end
    
    # Plot the curve with arrows to show increasing amplitude direction
    plot!(fig, re_vals, im_vals, label="-1/N(A): $label", lw=2.5, arrow=true, hover=A_range)
    
    # Mark the start (low A) and end (high A)
    scatter!(fig, [re_vals[1]], [im_vals[1]], label="Low A", markershape=:circle)
    scatter!(fig, [re_vals[end]], [im_vals[end]], label="High A", markershape=:square)
    
    return fig
end

# ==============================================================================
## Examples of usage
# ==============================================================================
# 1. Sign (Relay) function: f(x) = sgn(x)
relay(x) = sign(x)
println("Relay N(2): ", compute_describing_function(relay, 2.0)) 
# Theoretical check: 4/(π*A) = 4/(2π) ≈ 0.6366

# 2. Saturation function
saturation(x, limit=1.0) = clamp(x, -limit, limit)
println("Saturation N(2): ", compute_describing_function(x -> saturation(x, 1.0), 2.0))

# 3. Dead-zone function
deadzone(x, d=0.1) = abs(x) < d ? 0.0 : x - sign(x)*d
println("Dead-zone N(2): ", compute_describing_function(x -> deadzone(x, 0.1), 2.0))


# ==============================================================================
## Nyquist Plot with Describing Function Overlay
# ==============================================================================

using ControlSystems, Plots

nonlin = saturation

s = tf("s")
G = 10 / (s^3 + 2s^2 + s + 1)
# Create the linear Nyquist plot first
p_nyq = nyquistplot(G)
# Overlay the describing function
plot_negative_inverse_df(nonlin, label="deadzone"; fig=p_nyq,
    A_range=0.01:0.01:100)
# The intersection between the Nyquist plot and -1/N(A) indicates potential limit cycles. The amplitude of the limit cycle can be estimated from the corresponding A value, use the plotly() backend to display the amplitude on hover.

# ==============================================================================
## Simulation
# We can simulate the closed-loop system with the nonlinearity to verify the limit cycle prediction.
# ==============================================================================

nl = nonlinearity(nonlin)
sys_cl = feedback(G*nl)
plot(step(sys_cl, 200))