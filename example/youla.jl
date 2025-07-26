# 1. SETUP: Load packages and define the system
using ControlSystemsBase, JuMP, Ipopt, Plots, Printf

# Define the continuous-time plant model (a simple damped oscillator)
P_s = tf(10, [1, 2, 10])

# Define performance weights
# Weight for sensitivity S (good tracking at low frequencies)
Ws = tf(1, [1, 1e-4]) 
# Weight for control effort U = Q*S (penalize high-frequency control)
Wu = tf([1, 0], [1/20, 1])

# T = PQ
# CS = Q


# Discretize the system with a sample time Ts
Ts = 0.05
P = c2d(P_s, Ts)
Ws_d = c2d(Ws, Ts)
Wu_d = c2d(Wu, Ts)

# ---
# 2. YOULA PARAMETER (Q) SETUP
# We parameterize Q as an FIR filter of order Nq.
# Q(z) = q_0 + q_1*z⁻¹ + ... + q_{Nq-1}*z⁻⁽ᴺ𝑞⁻¹⁾
# The coefficients q_i are our optimization variables.
Nq = 15 # Order of the Q filter

# ---
# 3. OPTIMIZATION PROBLEM SETUP
# Frequency grid for optimization
w = [
    1e-12; # To get integral action (approximately)
    exp10.(LinRange(-2, log10(pi/Ts), 200))
]

# Evaluate transfer functions at each frequency in the grid
P_fr = freqrespv(P, w)
Ws_fr = freqrespv(Ws_d, w)
Wu_fr = freqrespv(Wu_d, w)

# Define the optimization model
model = Model(Ipopt.Optimizer)
set_silent(model) # Suppress solver output

# Define the decision variables (the FIR coefficients of Q)
@variable(model, q[1:Nq])

# Build expressions for Q(jω) in terms of the variables q
# Q(e^{jωT}) = Σ q_k * e^{-jωT(k-1)}
# We need to handle complex numbers by separating real and imaginary parts.
z = exp.(-im .* w .* Ts .* (0:Nq-1)') # [Nω x Nq] matrix
Qr = real(z) * q # Real part of Q at each frequency ω
Qi = imag(z) * q # Imaginary part of Q at each frequency ω

# Objective Function: Minimize ∫ |Ws*S|² + |Wu*Q|² dω
# which we approximate as a sum over our frequency grid.
# S = 1 - P*Q ✔
Sr = 1 .- (real.(P_fr) .* Qr .- imag.(P_fr) .* Qi)
Si = 0 .- (real.(P_fr) .* Qi .+ imag.(P_fr) .* Qr)

# Weighted sensitivity term |Ws*S|²
J_S = @expression(model, abs2.(Ws_fr) .* (Sr.^2 .+ Si.^2))

# Weighted control effort term |Wu*Q|²
J_U = @expression(model, abs2.(Wu_fr) .* (Qr.^2 .+ Qi.^2))

# Set the nonlinear objective function (sum over all frequencies)
@objective(model, Min, sum(J_S[i] + J_U[i] for i in eachindex(w)))

# Constraint: Peak of sensitivity function |S(jω)| ≤ γ
gamma_S = 1.5 # Maximum allowed sensitivity peak
@constraint(model, con[i=1:length(w)], Sr[i]^2 + Si[i]^2 <= gamma_S^2)

# ---
# 4. SOLVE AND RECONSTRUCT
println("Optimizing controller...")
optimize!(model)

# Check solver status
if termination_status(model) == MOI.LOCALLY_SOLVED
    println("Optimization successful.")
    # Extract optimal Q coefficients
    q_opt = value.(q)

    # Reconstruct the optimal Q and the final controller C
    Q_opt = tf(q_opt, [1; zeros(Nq-1)], Ts) |> ss
    C = minreal(Q_opt / (1 - P*Q_opt)) # C = Q / (1 - P*Q) ✔
else
    error("Optimization failed with status: $(termination_status(model))")
end

println("\nFinal Controller C(z):")
display(C)

# ---
# 5. ANALYSIS AND VISUALIZATION
println("\nAnalyzing closed-loop performance...")
S = sensitivity(P, C)     # Sensitivity
T = comp_sensitivity(P, C)   # Complementary Sensitivity

# Plot frequency responses
bodeplot(
    [S, T],
    w[2:end],
    lab=["Sensitivity S" "Comp. Sensitivity T"],
    title="Closed-Loop Frequency Response",
    linewidth=2
)
# Add the constraint line to the plot
hline!([gamma_S], lab="Constraint (γ = $(gamma_S))", color=:red, ls=:dash)

# Plot Gang of Four for a complete picture
gangoffourplot(P, C, title="Gang of Four", linewidth=2)

# Simulate the step response
# plot(step(T, 200), title="Closed-Loop Step Response", lab="y(t)", linewidth=2)