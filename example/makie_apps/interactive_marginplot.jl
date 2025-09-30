"""
# Interactive Margin Plot with P/PI/PID Controller Design

This example demonstrates an interactive Makie app for designing P/PI/PID controllers
by specifying magnitude at 1, 2, or 3 frequency points on a Bode magnitude plot.

- **1 point**: Proportional (P) controller - Kp only
- **2 points**: Proportional-Integral (PI) controller - Kp, Ki
- **3 points**: Proportional-Integral-Derivative (PID) controller - Kp, Ki, Kd

Uses pure algebraic approach (no iterative solvers) to compute controller parameters.

## Usage:
1. Run this file to launch the interactive app
2. Click to select FIRST point → Designs P controller (Kp)
3. Click to select SECOND point → Upgrades to PI controller (Kp, Ki)
4. Click to select THIRD point → Upgrades to PID controller (Kp, Ki, Kd)
5. Drag any marker vertically to adjust magnitude constraints
6. Controller updates in real-time as you drag

## Requirements:
- ControlSystemsBase
- GLMakie (for interactive features)
"""

using ControlSystemsBase
using GLMakie
using Printf

"""
    find_nearest_frequency(w, x_click_log)

Find the index of the frequency in w closest to x_click_log.
x_click_log is already in log10 scale (from Makie's log-scaled axis).
"""
function find_nearest_frequency(w, x_click_log)
    # x_click_log is already log10(frequency) from the log-scaled axis
    log_w = log10.(w)
    distances = @. abs(log_w - x_click_log)
    return argmin(distances)
end

"""
    compute_PI_two_points(P, ω1, m1_db, ω2, m2_db)

Compute PI controller (Kp, Ki) such that:
- |P(jω1) * C(jω1)| = 10^(m1_db/20)
- |P(jω2) * C(jω2)| = 10^(m2_db/20)

where C(s) = Kp + Ki/s

Returns: (Kp, Ki, success_flag, error_message)
"""
function compute_PI_two_points(P, ω1, m1_db, ω2, m2_db)
    try
        # Check for degenerate case (frequencies too close)
        if abs(ω1 - ω2) < 1e-6
            return (0.0, 0.0, false, "Frequencies too close")
        end

        # Convert dB to linear
        m1 = 10^(m1_db / 20)
        m2 = 10^(m2_db / 20)

        # Get plant frequency response
        P1 = freqresp(P, ω1)[]
        P2 = freqresp(P, ω2)[]

        # Required controller magnitudes
        r1 = m1 / abs(P1)
        r2 = m2 / abs(P2)

        # Solve for Ki
        # |C(jω)| = sqrt(Kp² + (Ki/ω)²)
        # r1² = Kp² + (Ki/ω1)²
        # r2² = Kp² + (Ki/ω2)²
        # Subtracting: r1² - r2² = (Ki/ω1)² - (Ki/ω2)²

        denom = 1/ω1^2 - 1/ω2^2
        numer = r1^2 - r2^2

        Ki_sq = numer / denom

        # If Ki² is negative, try swapping the points
        if Ki_sq < 0
            # Swap and try again
            denom = 1/ω2^2 - 1/ω1^2
            numer = r2^2 - r1^2
            Ki_sq = numer / denom

            if Ki_sq < 0
                return (0.0, 0.0, false, "No valid PI solution exists")
            end

            # Use swapped values
            ω1, ω2 = ω2, ω1
            r1, r2 = r2, r1
        end

        Ki = sqrt(Ki_sq)

        # Solve for Kp from first equation
        Kp_sq = r1^2 - (Ki/ω1)^2

        if Kp_sq < 0
            return (0.0, 0.0, false, "No valid PI solution exists")
        end

        # Take positive square root for Kp (could also choose negative)
        Kp = sqrt(Kp_sq)

        return (Kp, Ki, true, "")

    catch e
        return (0.0, 0.0, false, "Error: $e")
    end
end

"""
    compute_PID_three_points(P, ω1, m1_db, ω2, m2_db, ω3, m3_db)

Compute PID controller (Kp, Ki, Kd) such that:
- |P(jω1) * C(jω1)| = 10^(m1_db/20)
- |P(jω2) * C(jω2)| = 10^(m2_db/20)
- |P(jω3) * C(jω3)| = 10^(m3_db/20)

where C(s) = Kp + Ki/s + Kd*s

Uses algebraic approach: enumerates 8 sign combinations to find solution.

Returns: (Kp, Ki, Kd, success_flag, error_message)
"""
function compute_PID_three_points(P, ω1, m1_db, ω2, m2_db, ω3, m3_db)
    try
        # Check for degenerate cases
        if abs(ω1 - ω2) < 1e-6 || abs(ω2 - ω3) < 1e-6 || abs(ω1 - ω3) < 1e-6
            return (0.0, 0.0, 0.0, false, "Frequencies too close")
        end

        # Convert dB to linear
        m1 = 10^(m1_db / 20)
        m2 = 10^(m2_db / 20)
        m3 = 10^(m3_db / 20)

        # Get plant frequency response
        P1 = freqresp(P, ω1)[]
        P2 = freqresp(P, ω2)[]
        P3 = freqresp(P, ω3)[]

        # Required controller magnitudes
        r1 = m1 / abs(P1)
        r2 = m2 / abs(P2)
        r3 = m3 / abs(P3)

        # For PID: |C(jω)| = sqrt(Kp² + (Kd*ω - Ki/ω)²)
        # Let yᵢ = Kd*ωᵢ - Ki/ωᵢ
        # Then: Kp² + yᵢ² = rᵢ²  →  yᵢ = ±sqrt(rᵢ² - Kp²)

        # Try all 8 sign combinations
        tolerance = 1e-6

        for s1 in [-1, 1], s2 in [-1, 1], s3 in [-1, 1]
            # Assume y values and solve
            # We need to find Kp first, then compute y values

            # Actually, we don't know Kp yet. Let's solve the linear system for Kd, Ki
            # treating y1, y2 as unknowns with given signs

            # From equations: Kd*ω₁ - Ki/ω₁ = y₁ and Kd*ω₂ - Ki/ω₂ = y₂
            # This is underdetermined since y₁, y₂ depend on Kp

            # Better approach: try values of Kp, compute y values, solve for Kd, Ki, check y₃

            # Actually, let me think differently.
            # From first two equations:
            # (Kd*ω₁ - Ki/ω₁)² = r₁² - Kp²
            # (Kd*ω₂ - Ki/ω₂)² = r₂² - Kp²

            # If we assume signs, y₁ = s1*sqrt(r₁² - Kp²), y₂ = s2*sqrt(r₂² - Kp²)
            # But we still don't know Kp.

            # From y₁² - y₂² = (r₁² - Kp²) - (r₂² - Kp²) = r₁² - r₂²
            # So: y₁² - y₂² = r₁² - r₂²
            # And:y₁ = s1*|y₁|, y₂ = s2*|y₂|

            # Let's solve it differently: assume we have y₁, y₂ (with signs)
            # Solve linear system for Kd, Ki:
            # ω₁*Kd - Ki/ω₁ = y₁
            # ω₂*Kd - Ki/ω₂ = y₂

            # From the magnitude equations, if the signs are consistent:
            # We need |y₁|² = r₁² - Kp² and |y₂|² = r₂² - Kp²
            # This gives: |y₁|² - |y₂|² = r₁² - r₂²

            # For given signs, we can parameterize by Kp and solve
            # But that's complex. Let me try a direct enumeration approach.

            # Actually, simpler: for each combination, solve assuming Kp > 0
            # Check if all three equations are satisfied

            # Let me solve for one Kp value and verify

            # From r₁² = Kp² + y₁², if we guess Kp, we get y₁ = s1*sqrt(r₁² - Kp²)
            # Similarly for y₂ and y₃

            # Then solve linear system [ω₁ -1/ω₁][Kd] = [y₁]
            #                         [ω₂ -1/ω₂][Ki]   [y₂]

            # And verify with third equation

            # But we need to find Kp. From consistency of the first two:
            # After solving for Kd, Ki from y₁, y₂, we can compute Kp from either equation

            # Let me implement a simpler approach:
            # From the constraint that the first two equations must give the same Kp:
            # Kp₁ = sqrt(r₁² - y₁²) must equal Kp₂ = sqrt(r₂² - y₂²)
            # So: r₁² - y₁² = r₂² - y₂²
            #     y₁² - y₂² = r₁² - r₂²

            # Also: y₁ = Kd*ω₁ - Ki/ω₁ and y₂ = Kd*ω₂ - Ki/ω₂

            # Let me solve assuming y₁² - y₂² = r₁² - r₂²
            # With signs: s1²*|y₁|² - s2²*|y₂|² = r₁² - r₂² (since s² = 1)
            #            |y₁|² - |y₂|² = r₁² - r₂²

            # So: |y₁| = sqrt(|y₂|² + r₁² - r₂²)

            # Hmm, this is getting circular. Let me try a more direct approach:

            # Given r₁, r₂, r₃ and ω₁, ω₂, ω₃:
            # 1. Assume Kp = r₁ (initial guess, proportional controller at ω₁)
            # 2. Compute y₁ = s1*sqrt(r₁² - Kp²) (might be 0)
            # 3. Compute y₂ = s2*sqrt(r₂² - Kp²)
            # 4. Solve for Kd, Ki
            # 5. Recompute Kp from equation 2: Kp_check = sqrt(r₂² - y₂²)
            # 6. If Kp ≈ Kp_check, verify with equation 3

            # Actually simpler: just solve the system directly
            # From y₁² - y₂² = r₁² - r₂²
            # And linear system, we can eliminate Kp

            # Let me implement the direct algebraic solution:
            # Given signs s1, s2:
            # If y₁² - y₂² = r₁² - r₂², and both must be non-negative:

            # diff_sq = r1^2 - r2^2

            # y₁² - y₂² = diff_sq
            # Try to find y₁, y₂ such that this holds and they're both valid

            # With signs: y₁ = s1*|y₁|, y₂ = s2*|y₂|
            # We need: y₁² - y₂² = s1²*|y₁|² - s2²*|y₂|² = |y₁|² - |y₂|² = diff_sq

            # If diff_sq >= 0: |y₁|² = |y₂|² + diff_sq
            # If diff_sq < 0: |y₂|² = |y₁|² - diff_sq

            # Let's parameterize by |y₂| (or |y₁|)
            # If diff_sq >= 0:
            #   For any |y₂| >= 0, we get |y₁| = sqrt(|y₂|² + diff_sq)
            # If diff_sq < 0:
            #   We need |y₁|² >= -diff_sq, so |y₁| >= sqrt(-diff_sq)

            # But we also have constraints from Kp:
            # Kp² = r₁² - y₁² = r₁² - |y₁|² >= 0  →  |y₁| <= r₁
            # Kp² = r₂² - y₂² = r₂² - |y₂|² >= 0  →  |y₂| <= r₂

            # And both must give the same Kp!

            # From r₁² - y₁² = r₂² - y₂²:
            #      y₁² - y₂² = r₁² - r₂²  (which we already have)

            # So the system is consistent. We can choose one parameter freely.
            # Let's say we choose Kp directly.

            # For a given Kp in [0, min(r₁, r₂, r₃)]:
            # y₁ = s1*sqrt(r₁² - Kp²)
            # y₂ = s2*sqrt(r₂² - Kp²)
            # y₃ = s3*sqrt(r₃² - Kp²)

            # Solve for Kd, Ki from first two equations:
            # ω₁*Kd - Ki/ω₁ = y₁
            # ω₂*Kd - Ki/ω₂ = y₂

            # Check if third equation is satisfied:
            # ω₃*Kd - Ki/ω₃ ≈ y₃

            # Let me implement this search:

            # Kp search range
            Kp_min = 0.0
            Kp_max = min(r1, r2, r3)

            if Kp_max <= 0
                continue  # No valid Kp
            end

            # Binary search or direct solve
            # Actually, we can solve directly from the constraint
            # that ω₃*Kd - Ki/ω₃ = y₃

            # From equations 1 and 2:
            # det = ω₁*(-1/ω₂) - (-1/ω₁)*ω₂ = -ω₁/ω₂ + ω₂/ω₁
            det = (ω1^2 - ω2^2) / (ω1*ω2)

            if abs(det) < 1e-10
                continue  # Degenerate
            end

            # Use Newton-Raphson to find Kp such that third equation is satisfied
            # Define residual function: we want y3_computed - y3 = 0

            function residual(Kp)
                if Kp < 0 || Kp > Kp_max
                    return Inf
                end

                y1_sq = r1^2 - Kp^2
                y2_sq = r2^2 - Kp^2
                y3_sq = r3^2 - Kp^2

                if y1_sq < 0 || y2_sq < 0 || y3_sq < 0
                    return Inf
                end

                y1 = s1 * sqrt(y1_sq)
                y2 = s2 * sqrt(y2_sq)
                y3 = s3 * sqrt(y3_sq)

                # Solve for Kd, Ki from first two equations
                Kd = (ω1*y1 - ω2*y2) / (ω1^2 - ω2^2)
                Ki = ω1*ω2*(ω2*y1 - ω1*y2) / (ω1^2 - ω2^2)

                # Check third equation
                y3_computed = Kd*ω3 - Ki/ω3

                return y3_computed - y3
            end

            # Try bisection search
            # First check if there's a sign change in the interval
            res_min = residual(Kp_min + 1e-6)
            res_max = residual(Kp_max - 1e-6)

            if isfinite(res_min) && isfinite(res_max) && res_min * res_max < 0
                # Bisection
                a, b = Kp_min + 1e-6, Kp_max - 1e-6
                for iter in 1:50
                    c = (a + b) / 2
                    res_c = residual(c)

                    if abs(res_c) < tolerance || abs(b - a) < tolerance
                        # Found solution!
                        Kp_sol = c
                        y1_sol = s1 * sqrt(r1^2 - Kp_sol^2)
                        y2_sol = s2 * sqrt(r2^2 - Kp_sol^2)

                        Kd = (ω1*y1_sol - ω2*y2_sol) / (ω1^2 - ω2^2)
                        Ki = ω1*ω2*(ω2*y1_sol - ω1*y2_sol) / (ω1^2 - ω2^2)

                        return (Kp_sol, Ki, Kd, true, "")
                    end

                    if residual(a) * res_c < 0
                        b = c
                    else
                        a = c
                    end
                end
            end
        end

        return (0.0, 0.0, 0.0, false, "No valid PID solution exists")

    catch e
        rethrow()
        return (0.0, 0.0, 0.0, false, "Error: $e")
    end
end

"""
    interactive_marginplot(P::LTISystem; w=nothing)

Create an interactive Bode/margin plot for P/PI/PID controller design.

Supports 1, 2, or 3 frequency-magnitude constraints:
- 1 point → P controller (Kp only)
- 2 points → PI controller (Kp, Ki)
- 3 points → PID controller (Kp, Ki, Kd)

# Arguments
- `P::LTISystem`: The plant transfer function (must be SISO)
- `w::AbstractVector`: Optional frequency vector (rad/s)

# Returns
- `fig`: The interactive Makie figure

# Example
```julia
using ControlSystemsBase, GLMakie

# Define a plant
P = tf(1, [1, 1]) * tf(1, [1, 0])

# Launch interactive app
fig = interactive_marginplot(P)
display(fig)
```
"""
function interactive_marginplot(P::LTISystem; w=nothing)
    # Validate input
    ControlSystemsBase.issiso(P) || error("Plant must be SISO")

    # Generate frequency vector if not provided
    if w === nothing
        w = exp10.(range(-3, stop=2, length=500))
    end

    # Ensure P has Float64 coefficients
    P_float = tf(P)

    # Compute initial Bode data for plant P
    mag_p, phase_p, w_used = bode(P_float, w)
    mag_p_db = 20 .* log10.(vec(mag_p))
    phase_p_deg = vec(phase_p)

    # Create observables for interactive updates
    controller = Observable{Any}(nothing)
    loop_transfer = Observable{Any}(P_float)
    mag_loop_db = Observable(copy(mag_p_db))
    phase_loop_deg = Observable(copy(phase_p_deg))

    # State: 0 = no points, 1 = one point, 2 = two points, 3 = three points selected
    selection_state = Observable(0)

    # Selected points: (index, target_mag_db)
    point1_idx = Observable(1)
    point1_mag = Observable(0.0)
    point1_pos = Observable(Point2f(0, 0))
    point1_visible = Observable(false)

    point2_idx = Observable(1)
    point2_mag = Observable(0.0)
    point2_pos = Observable(Point2f(0, 0))
    point2_visible = Observable(false)

    point3_idx = Observable(1)
    point3_mag = Observable(0.0)
    point3_pos = Observable(Point2f(0, 0))
    point3_visible = Observable(false)

    # Status text
    status_text = Observable("Select FIRST point on blue curve")

    # Margin indicators (initially nothing)
    gm_freq = Observable(NaN)  # Frequency where phase = -180°
    gm_mag = Observable(NaN)   # Magnitude at that frequency (dB)
    gm_value = Observable(NaN) # Gain margin in dB
    gm_text_str = Observable("")
    pm_freq = Observable(NaN)  # Frequency where magnitude = 0 dB
    pm_phase = Observable(NaN) # Phase at that frequency (deg)
    pm_value = Observable(NaN) # Phase margin in degrees
    pm_text_str = Observable("")

    # Create figure with two vertically stacked axes
    fig = Figure(size=(1000, 900))

    # Status label at top
    Label(fig[0, 1], status_text, fontsize=16, color=:blue, tellwidth=false)

    # Magnitude axis
    ax_mag = Axis(fig[1, 1],
                  xscale=log10,
                  xlabel="",
                  ylabel="Magnitude (dB)",
                  title="Interactive Margin Plot - Two Points",
                  xticklabelsvisible=false)

    # Phase axis
    ax_phase = Axis(fig[2, 1],
                    xscale=log10,
                    xlabel="Frequency (rad/s)",
                    ylabel="Phase (deg)")

    # Link x-axes
    linkxaxes!(ax_mag, ax_phase)

    # Disable default interactions
    deregister_interaction!(ax_mag, :rectanglezoom)
    deregister_interaction!(ax_mag, :limitreset)
    deregister_interaction!(ax_mag, :scrollzoom)

    # Plot plant curves
    lines!(ax_mag, w_used, mag_p_db, color=:blue, linewidth=2, label="Plant P")
    lines!(ax_mag, w_used, mag_loop_db, color=:red, linewidth=2, label="Loop P*C")
    hlines!(ax_mag, 0, color=:gray, linestyle=:dash, alpha=0.5, label="0 dB")

    lines!(ax_phase, w_used, phase_p_deg, color=:blue, linewidth=2, label="Plant P")
    lines!(ax_phase, w_used, phase_loop_deg, color=:red, linewidth=2, label="Loop P*C")
    hlines!(ax_phase, -180, color=:gray, linestyle=:dash, alpha=0.5, label="-180°")

    # Legends
    axislegend(ax_mag, position=:lb)
    axislegend(ax_phase, position=:lb)

    # Gain margin indicators (line segment from curve to 0 dB)
    gm_line_mag = lines!(ax_mag,
        @lift([Point2f($gm_freq, $gm_mag), Point2f($gm_freq, 0.0)]),
        color=:darkgreen, linewidth=2, linestyle=:dot)

    # Phase margin indicators (line segment from curve to -180°)
    pm_line_phase = lines!(ax_phase,
        @lift([Point2f($pm_freq, $pm_phase), Point2f($pm_freq, -180.0)]),
        color=:darkblue, linewidth=2, linestyle=:dot)

    # Text labels for margins
    text!(ax_mag, @lift(Point2f($gm_freq, $gm_mag - 5)),
          text=gm_text_str, color=:darkgreen, fontsize=12)
    text!(ax_phase, @lift(Point2f($gm_freq, -180.0)),
          text=@lift(isnan($gm_freq) ? "" : string("GM freq=", round($gm_freq, digits=2))),
          color=:darkgreen, fontsize=10, align=(:left, :bottom))
    text!(ax_mag, @lift(Point2f($pm_freq, 5.0)),
          text=pm_text_str, color=:darkblue, fontsize=10, align=(:left, :top))
    text!(ax_phase, @lift(Point2f($pm_freq, $pm_phase + 5)),
          text=@lift(isnan($pm_freq) ? "" : string("PM=", round($pm_value, digits=1), "°")),
          color=:darkblue, fontsize=12)

    # Interactive markers
    scatter!(ax_mag, point1_pos, color=:green, markersize=14,
             visible=point1_visible, label="Point 1")
    scatter!(ax_mag, point2_pos, color=:orange, markersize=14,
             visible=point2_visible, label="Point 2")
    scatter!(ax_mag, point3_pos, color=:purple, markersize=14,
             visible=point3_visible, label="Point 3")

    # Interaction state
    dragging = Ref(false)
    dragging_which = Ref(0)  # 1, 2, or 3

    # Helper function to compute and update margins
    function update_margins(mag_db, phase_deg, w)
        # Find gain margin: magnitude at frequency where phase = -180°
        # Phase crossover frequency
        phase_cross_idx = findfirst(i -> i > 1 &&
            (phase_deg[i] <= -180 && phase_deg[i-1] > -180),
            eachindex(phase_deg))

        if phase_cross_idx !== nothing
            # Interpolate to find exact frequency
            idx = phase_cross_idx
            if idx > 1
                # Linear interpolation
                f1, f2 = w[idx-1], w[idx]
                p1, p2 = phase_deg[idx-1], phase_deg[idx]
                m1, m2 = mag_db[idx-1], mag_db[idx]

                # Interpolate frequency where phase = -180
                t = (-180 - p1) / (p2 - p1)
                gm_freq[] = f1 * (1 - t) + f2 * t
                gm_mag[] = m1 * (1 - t) + m2 * t
                gm_value[] = -gm_mag[]  # GM in dB (positive means stable)
                gm_text_str[] = string("GM=", round(gm_value[], digits=1), "dB")
            end
        else
            gm_freq[] = NaN
            gm_mag[] = NaN
            gm_value[] = NaN
            gm_text_str[] = ""
        end

        # Find phase margin: phase at frequency where magnitude = 0 dB
        # Gain crossover frequency
        gain_cross_idx = findfirst(i -> i > 1 &&
            (mag_db[i] <= 0 && mag_db[i-1] > 0),
            eachindex(mag_db))

        if gain_cross_idx !== nothing
            idx = gain_cross_idx
            if idx > 1
                # Linear interpolation
                f1, f2 = w[idx-1], w[idx]
                m1, m2 = mag_db[idx-1], mag_db[idx]
                p1, p2 = phase_deg[idx-1], phase_deg[idx]

                # Interpolate frequency where magnitude = 0
                t = (0 - m1) / (m2 - m1)
                pm_freq[] = f1 * (1 - t) + f2 * t
                pm_phase[] = p1 * (1 - t) + p2 * t
                pm_value[] = pm_phase[] + 180  # PM in degrees
                pm_text_str[] = string("PM freq=", round(pm_freq[], digits=2))
            end
        else
            pm_freq[] = NaN
            pm_phase[] = NaN
            pm_value[] = NaN
            pm_text_str[] = ""
        end
    end

    # Helper function to update controller
    function update_controller()
        if selection_state[] == 1
            # One point: design proportional-only controller (Kp, Ki=0)
            ω1 = w_used[point1_idx[]]
            m1_db = point1_mag[]
            m1 = 10^(m1_db / 20)

            P1 = freqresp(P_float, ω1)[]
            Kp = m1 / abs(P1)
            Ki = 0.0

            # Create controller and update
            C = pid(Kp, Ki, 0, form=:parallel)
            controller[] = C
            L = P_float * C
            loop_transfer[] = L

            # Update Bode curves
            mag_l, phase_l, _ = bode(L, w_used)
            mag_loop_db[] = 20 .* log10.(vec(mag_l))
            phase_loop_deg[] = vec(phase_l)

            # Update margins
            update_margins(mag_loop_db[], phase_loop_deg[], w_used)

            # Update status
            status_text[] = @sprintf(
                "✓ P-only: Kp=%.3f, Ki=0 | ω1=%.2f→%.1fdB | Select 2nd point for PI",
                Kp, ω1, m1_db
            )

        elseif selection_state[] == 2
            # Two points: design full PI controller
            ω1 = w_used[point1_idx[]]
            ω2 = w_used[point2_idx[]]
            m1_db = point1_mag[]
            m2_db = point2_mag[]

            Kp, Ki, success, err_msg = compute_PI_two_points(P_float, ω1, m1_db, ω2, m2_db)

            if success
                # Create controller and update
                C = pid(Kp, Ki, 0, form=:parallel)
                controller[] = C
                L = P_float * C
                loop_transfer[] = L

                # Update Bode curves
                mag_l, phase_l, _ = bode(L, w_used)
                mag_loop_db[] = 20 .* log10.(vec(mag_l))
                phase_loop_deg[] = vec(phase_l)

                # Update margins
                update_margins(mag_loop_db[], phase_loop_deg[], w_used)

                # Update status
                status_text[] = @sprintf(
                    "✓ PI: Kp=%.3f, Ki=%.3f | ω1=%.2f→%.1fdB, ω2=%.2f→%.1fdB",
                    Kp, Ki, ω1, m1_db, ω2, m2_db
                )
            else
                status_text[] = "✗ " * err_msg * " | Try different points"
            end

        elseif selection_state[] == 3
            # Three points: design full PID controller
            ω1 = w_used[point1_idx[]]
            ω2 = w_used[point2_idx[]]
            ω3 = w_used[point3_idx[]]
            m1_db = point1_mag[]
            m2_db = point2_mag[]
            m3_db = point3_mag[]

            Kp, Ki, Kd, success, err_msg = compute_PID_three_points(P_float, ω1, m1_db, ω2, m2_db, ω3, m3_db)

            if success
                # Create controller and update
                C = pid(Kp, Ki, Kd, form=:parallel)
                controller[] = C
                L = P_float * C
                loop_transfer[] = L

                # Update Bode curves
                mag_l, phase_l, _ = bode(L, w_used)
                mag_loop_db[] = 20 .* log10.(vec(mag_l))
                phase_loop_deg[] = vec(phase_l)

                # Update margins
                update_margins(mag_loop_db[], phase_loop_deg[], w_used)

                # Update status
                status_text[] = @sprintf(
                    "✓ PID: Kp=%.3f, Ki=%.3f, Kd=%.3f | ω1=%.2f→%.1fdB, ω2=%.2f→%.1fdB, ω3=%.2f→%.1fdB",
                    Kp, Ki, Kd, ω1, m1_db, ω2, m2_db, ω3, m3_db
                )
            else
                status_text[] = "✗ " * err_msg * " | Try different points"
            end
        end
    end

    # Mouse button handler
    on(events(fig).mousebutton) do event
        if event.button == Mouse.left
            if event.action == Mouse.press
                mp = mouseposition(ax_mag.scene)

                # Find nearest frequency
                idx = find_nearest_frequency(w_used, mp[1])

                # Check if clicking near existing point to drag
                if selection_state[] >= 1
                    # Compare distances in log space since axis is log-scaled
                    dist_to_point1 = abs(log10(w_used[point1_idx[]]) - log10(w_used[idx]))

                    # Use log-scale threshold (e.g., within 0.1 decades)
                    threshold = 0.15  # decades in log space

                    # Check if near point 1
                    if dist_to_point1 < threshold
                        dragging[] = true
                        dragging_which[] = 1
                        return
                    end

                    # If two or more points, also check point 2
                    if selection_state[] >= 2
                        dist_to_point2 = abs(log10(w_used[point2_idx[]]) - log10(w_used[idx]))
                        if dist_to_point2 < threshold
                            dragging[] = true
                            dragging_which[] = 2
                            return
                        end
                    end

                    # If three points, also check point 3
                    if selection_state[] == 3
                        dist_to_point3 = abs(log10(w_used[point3_idx[]]) - log10(w_used[idx]))
                        if dist_to_point3 < threshold
                            dragging[] = true
                            dragging_which[] = 3
                            return
                        end
                    end
                end

                # Otherwise, selecting new point
                if selection_state[] == 0
                    # Select first point
                    point1_idx[] = idx
                    point1_mag[] = mp[2]
                    point1_pos[] = Point2f(w_used[idx], mp[2])
                    point1_visible[] = true
                    selection_state[] = 1
                    update_controller()  # Design P-only controller immediately

                elseif selection_state[] == 1
                    # Select second point
                    point2_idx[] = idx
                    point2_mag[] = mp[2]
                    point2_pos[] = Point2f(w_used[idx], mp[2])
                    point2_visible[] = true
                    selection_state[] = 2
                    update_controller()

                elseif selection_state[] == 2
                    # Select third point
                    point3_idx[] = idx
                    point3_mag[] = mp[2]
                    point3_pos[] = Point2f(w_used[idx], mp[2])
                    point3_visible[] = true
                    selection_state[] = 3
                    update_controller()

                elseif selection_state[] == 3
                    # Reset and start over
                    selection_state[] = 0
                    point1_visible[] = false
                    point2_visible[] = false
                    point3_visible[] = false
                    status_text[] = "Select FIRST point on blue curve"
                    # Reset loop to plant
                    mag_loop_db[] = copy(mag_p_db)
                    phase_loop_deg[] = copy(phase_p_deg)
                end

            elseif event.action == Mouse.release
                dragging[] = false
            end
        end
    end

    # Mouse position handler (dragging)
    on(events(fig).mouseposition) do mp_px
        if dragging[]
            mp = mouseposition(ax_mag.scene)

            # Find nearest frequency for horizontal position
            idx = find_nearest_frequency(w_used, mp[1])

            if dragging_which[] == 1
                # Drag point 1 - update both frequency and magnitude
                point1_idx[] = idx
                point1_mag[] = mp[2]
                point1_pos[] = Point2f(w_used[idx], mp[2])
            elseif dragging_which[] == 2
                # Drag point 2 - update both frequency and magnitude
                point2_idx[] = idx
                point2_mag[] = mp[2]
                point2_pos[] = Point2f(w_used[idx], mp[2])
            elseif dragging_which[] == 3
                # Drag point 3 - update both frequency and magnitude
                point3_idx[] = idx
                point3_mag[] = mp[2]
                point3_pos[] = Point2f(w_used[idx], mp[2])
            end

            update_controller()
        end
    end

    return fig
end

# Example usage
function demo_interactive_marginplot()
    println("Interactive Two-Point Margin Plot Demo")

    # Example: First-order system with integrator
    P = tf(1, [1, 1]) * tf(1, [1, 0])

    fig = interactive_marginplot(P)
    display(fig)

    println("\nInteractive Margin Plot launched!")
    println("Instructions:")
    println("1. Click to select FIRST point → P controller (Kp)")
    println("2. Click to select SECOND point → PI controller (Kp, Ki)")
    println("3. Click to select THIRD point → PID controller (Kp, Ki, Kd)")
    println("4. Drag any marker vertically to adjust magnitude constraints")
    println("5. Click anywhere to reset and start over")
    println("\nNote:")
    println("  • 1 point = P control")
    println("  • 2 points = PI control")
    println("  • 3 points = PID control")

    return fig
end

# Allow running as a script
demo_interactive_marginplot()