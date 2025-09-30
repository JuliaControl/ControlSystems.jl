"""
# Interactive Margin Plot with Two-Point PI Controller Design

This example demonstrates an interactive Makie app for designing PI controllers
by specifying magnitude at TWO frequency points on a Bode magnitude plot.
This properly utilizes both degrees of freedom (Kp and Ki) of a PI controller.

## Usage:
1. Run this file to launch the interactive app
2. Click to select FIRST point → Designs proportional controller (Kp only, Ki=0)
3. Drag green marker vertically to adjust gain
4. Click to select SECOND point → Upgrades to full PI controller (Kp + Ki)
5. Drag either marker to adjust constraints, controller updates in real-time

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
    interactive_marginplot(P::LTISystem; w=nothing)

Create an interactive Bode/margin plot for two-point PI controller design.

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

    # State: 0 = no points, 1 = one point, 2 = two points selected
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

    # Status text
    status_text = Observable("Select FIRST point on blue curve")

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

    # Interactive markers
    scatter!(ax_mag, point1_pos, color=:green, markersize=14,
             visible=point1_visible, label="Point 1")
    scatter!(ax_mag, point2_pos, color=:orange, markersize=14,
             visible=point2_visible, label="Point 2")

    # Interaction state
    dragging = Ref(false)
    dragging_which = Ref(0)  # 1 or 2

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

                # Update status
                status_text[] = @sprintf(
                    "✓ PI: Kp=%.3f, Ki=%.3f | ω1=%.2f→%.1fdB, ω2=%.2f→%.1fdB",
                    Kp, Ki, ω1, m1_db, ω2, m2_db
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

                    # If two points, also check point 2
                    if selection_state[] == 2
                        dist_to_point2 = abs(log10(w_used[point2_idx[]]) - log10(w_used[idx]))
                        if dist_to_point2 < threshold
                            dragging[] = true
                            dragging_which[] = 2
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
                    # Reset and start over
                    selection_state[] = 0
                    point1_visible[] = false
                    point2_visible[] = false
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

            if dragging_which[] == 1
                # Drag point 1
                point1_mag[] = mp[2]
                point1_pos[] = Point2f(w_used[point1_idx[]], mp[2])
            elseif dragging_which[] == 2
                # Drag point 2
                point2_mag[] = mp[2]
                point2_pos[] = Point2f(w_used[point2_idx[]], mp[2])
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
    println("1. Click to select FIRST point → Proportional controller (Kp only)")
    println("2. Drag green marker to adjust gain")
    println("3. Click to select SECOND point → Full PI controller (Kp + Ki)")
    println("4. Drag either marker vertically to adjust magnitude")
    println("5. Click anywhere to reset and start over")
    println("\nNote: One point = P control, Two points = full PI control!")

    return fig
end

# Allow running as a script
demo_interactive_marginplot()