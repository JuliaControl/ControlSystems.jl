"""
# Interactive Nyquist Loop-Shaping with PI Controller Design

This example demonstrates an interactive Makie app for designing PI controllers
by dragging points on the Nyquist plot. When you grab and move a point on the
Nyquist curve, the app uses `loopshapingPI` to design a controller that moves
the loop transfer function through the desired location.

## Usage:
1. Run this file to launch the interactive app
2. Click near any point on the blue Nyquist curve (plant P)
3. Drag the point to where you want the loop transfer function (P*C) to be
4. Release to see the designed PI controller and updated red curve (P*C)
5. Controller parameters are displayed in the title

## Requirements:
- ControlSystemsBase
- GLMakie (for interactive features)
"""

using ControlSystemsBase
using GLMakie
using Printf

"""
    find_nearest_point(re_data, im_data, x, y)

Find the index of the point in (re_data, im_data) closest to (x, y).
"""
function find_nearest_point(re_data, im_data, x, y)
    distances = @. sqrt((re_data - x)^2 + (im_data - y)^2)
    return argmin(distances)
end

"""
    interactive_nyquist_loopshaping(P::LTISystem; w=nothing)

Create an interactive Nyquist plot for loop-shaping PI controller design.

# Arguments
- `P::LTISystem`: The plant transfer function (must be SISO)
- `w::AbstractVector`: Optional frequency vector (rad/s)

# Returns
- `fig`: The interactive Makie figure

# Example
```julia
using ControlSystemsBase, GLMakie

# Define a plant (e.g., double integrator with lag)
P = tf(1, [1, 1]) * tf(1, [1, 0])

# Launch interactive app
fig = interactive_nyquist_loopshaping(P)
display(fig)
```
"""
function interactive_nyquist_loopshaping(P::LTISystem; w=nothing, xlims=(-4, 1), ylims=(-4, 2))
    # Validate input
    ControlSystemsBase.issiso(P) || error("Plant must be SISO")

    # Generate frequency vector if not provided
    if w === nothing
        w = exp10.(range(-3, stop=2, length=500))
    end

    # Ensure P has Float64 coefficients to avoid type issues
    P_float = tf(P)  # Convert to ensure Float64 type

    # Compute initial Nyquist data for plant P
    re_p, im_p, w_used = nyquist(P_float, w)
    re_data = vec(re_p)
    im_data = vec(im_p)

    # Create observables for interactive updates (use Any to avoid type issues)
    controller = Observable{Any}(nothing)
    loop_transfer = Observable{Any}(P_float)
    re_loop = Observable(copy(re_data))
    im_loop = Observable(copy(im_data))
    controller_text = Observable("No controller designed yet\nClick and drag a point on the blue curve")

    # Create figure and axis
    fig = Figure(size=(900, 700))
    ax = Axis(fig[1, 1],
              aspect=DataAspect(),
              xlabel="Real",
              ylabel="Imaginary",
              title="Interactive Nyquist Loop-Shaping\n" * controller_text[],
              limits=(xlims, ylims))

    # Disable default zoom/pan interactions so our drag works
    deregister_interaction!(ax, :rectanglezoom)
    deregister_interaction!(ax, :limitreset)
    deregister_interaction!(ax, :scrollzoom)

    # Plot plant Nyquist curve (blue) - keep reference to original for display
    lines!(ax, re_data, im_data, color=:blue, linewidth=2, label="Plant P")

    # Plot loop transfer function curve (red) - initially same as plant
    loop_line = lines!(ax, re_loop, im_loop, color=:red, linewidth=2, label="Loop P*C")

    # Draw connecting lines between corresponding frequency points
    # Downsample to avoid too many lines (every 10th point)
    line_indices = 1:10:length(re_data)
    for idx in line_indices
        # Create observables for line endpoints
        line_points = @lift [Point2f(re_data[idx], im_data[idx]),
                             Point2f($re_loop[idx], $im_loop[idx])]
        lines!(ax, line_points, color=(:gray, 0.3), linewidth=0.5)
    end

    # Add reference lines
    vlines!(ax, 0, color=:gray, alpha=0.3, linewidth=0.5)
    hlines!(ax, 0, color=:gray, alpha=0.3, linewidth=0.5)

    # Critical point
    scatter!(ax, [-1], [0], marker=:xcross, markersize=15, color=:red, label="Critical point")

    # Add unit circle
    θ = range(0, 2π, length=100)
    lines!(ax, cos.(θ), sin.(θ), color=:gray, linestyle=:dash, alpha=0.5)

    # Scatter points for interaction (initially hidden, shown on hover/drag)
    selected_point = Observable(Point2f(0, 0))
    selected_point_visible = Observable(false)
    scatter!(ax, selected_point, color=:green, markersize=12, visible=selected_point_visible)

    # Legend
    axislegend(ax, position=:lt)

    # Update title when controller text changes
    on(controller_text) do text
        ax.title = "Interactive Nyquist Loop-Shaping\n" * text
    end

    # Interaction state
    dragging = Ref(false)
    selected_idx = Ref(0)
    selected_freq = Ref(0.0)

    # Mouse button event handler
    on(events(fig).mousebutton) do event
        if event.button == Mouse.left
            if event.action == Mouse.press
                # Get mouse position in data coordinates
                mp = mouseposition(ax.scene)

                # Check if click is within axis limits
                if mp[1] >= xlims[1] && mp[1] <= xlims[2] &&
                   mp[2] >= ylims[1] && mp[2] <= ylims[2]

                    # Find nearest point on plant curve
                    idx = find_nearest_point(re_data, im_data, mp[1], mp[2])

                    # Check if close enough to start dragging (within 0.3 units)
                    dist = sqrt((re_data[idx] - mp[1])^2 + (im_data[idx] - mp[2])^2)
                    if dist < 0.3
                        dragging[] = true
                        selected_idx[] = idx
                        selected_freq[] = w_used[idx]
                        selected_point[] = Point2f(mp[1], mp[2])
                        selected_point_visible[] = true
                    end
                end

            elseif event.action == Mouse.release
                if dragging[]
                    # Stop dragging
                    dragging[] = false
                    selected_point_visible[] = false
                end
            end
        end
    end

    # Mouse position event handler (for dragging)
    on(events(fig).mouseposition) do mp_px
        if dragging[]
            # Convert pixel coordinates to data coordinates
            mp = mouseposition(ax.scene)

            # Update selected point position and design controller in real-time
            if mp[1] >= xlims[1] && mp[1] <= xlims[2] &&
               mp[2] >= ylims[1] && mp[2] <= ylims[2]
                selected_point[] = Point2f(mp[1], mp[2])

                # Design controller in real-time as mouse moves
                rl = sqrt(mp[1]^2 + mp[2]^2)
                ϕl = atan(mp[2], mp[1])
                ω = selected_freq[]

                try
                    # Design PI controller using loopshapingPI
                    result = loopshapingPI(P_float, ω; ϕl=ϕl, rl=rl, form=:parallel)
                    C = result.C
                    kp = result.kp
                    ki = result.ki

                    # Update observables
                    controller[] = C
                    L = P_float * C
                    loop_transfer[] = L

                    # Compute new Nyquist data for loop transfer function
                    re_l, im_l, _ = nyquist(L, w_used)
                    re_loop[] = vec(re_l)
                    im_loop[] = vec(im_l)

                    # Update info text
                    controller_text[] = @sprintf(
                        "ω = %.3f rad/s | Kp = %.3f | Ki = %.3f\nTarget: |L| = %.3f, ∠L = %.1f°",
                        ω, kp, ki, rl, rad2deg(ϕl)
                    )

                catch e
                    # Silently ignore errors during dragging to avoid spam
                    # The user can try a different position
                end
            end
        end
    end

    return fig
end

# Example usage with different plants
function demo_interactive_nyquist()
    # Example 1: First-order system with integrator
    println("Example 1: First-order system with integrator")
    P1 = tf(1, [1, 1]) * tf(1, [1, 0])
    fig1 = interactive_nyquist_loopshaping(P1)
    display(fig1)

    println("\nInteractive Nyquist app launched!")
    println("Instructions:")
    println("1. Click near any point on the BLUE curve (plant P)")
    println("2. Drag to where you want the RED curve (loop P*C) to pass")
    println("3. Release to design the PI controller")
    println("4. The controller parameters will be shown in the title")
    println("\nTip: Try dragging points to locations with phase margin ~180° + 30-60°")
    println("     (e.g., locations at angles around -150° to -120°)")

    return fig1
end

# Allow running as a script
demo_interactive_nyquist()