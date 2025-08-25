module ControlSystemsBaseMakieExt

using ControlSystemsBase
using ControlSystemsBase.CSMakie
using Makie
using LinearAlgebra
using Printf

# Import necessary functions from ControlSystemsBase
using ControlSystemsBase: downsample, _processfreqplot, _default_freq_vector,
                          _same_io_dims, _get_plotlabel, _to1series,
                          SimResult, StepInfo, RootLocusResult,
                          poles, tzeros, bode, nyquist, sigma, margin, 
                          sisomargin, relative_gain_array, rlocus,
                          input_names, output_names, state_names, system_name,
                          iscontinuous, isdiscrete, issiso, isrational,
                          integrator_excess, balance_statespace, LTISystem,
                          _PlotScale, _PlotScaleFunc, _PlotScaleStr, _span  # Use existing plot scale settings

# Helper function to get y-scale transform for Makie
function get_yscale_transform()
    ControlSystemsBase._PlotScaleFunc == :log10 ? Makie.log10 : identity
end

# ====== PZMap (Pole-Zero Map) ======
function CSMakie.pzmap(systems::Union{LTISystem, AbstractVector{<:LTISystem}}; 
                       hz=false, kwargs...)
    systems_vec = systems isa AbstractVector ? systems : [systems]
    
    fig = Figure()
    ax = Axis(fig[1,1], aspect = DataAspect(),
             title = "Pole-zero map",
             xlabel = hz ? "Real [Hz]" : "Real",
             ylabel = hz ? "Imag [Hz]" : "Imag")
    
    # Add grid at zero
    vlines!(ax, 0, color=:gray, alpha=0.3, linewidth=0.5)
    hlines!(ax, 0, color=:gray, alpha=0.3, linewidth=0.5)
    
    scale_factor = hz ? 1/(2π) : 1
    
    for (i, system) in enumerate(systems_vec)
        p = poles(system) .* scale_factor
        z = tzeros(system) .* scale_factor
        
        # Plot zeros as circles
        if !isempty(z)
            scatter!(ax, real(z), imag(z), 
                    marker=:circle, markersize=15, 
                    color=(:transparent, 0.5), strokewidth=2, 
                    strokecolor=Cycled(i), label="Zeros")
        end
        
        # Plot poles as x's
        if !isempty(p)
            scatter!(ax, real(p), imag(p), 
                    marker=:xcross, markersize=15,
                    color=Cycled(i), label="Poles")
        end
        
        # Add unit circle for discrete systems
        if isdiscrete(system) && i == 1  # Only draw once
            θ = range(0, 2π, length=100)
            circle_scale = scale_factor
            lines!(ax, circle_scale .* cos.(θ), circle_scale .* sin.(θ), 
                  color=:gray, linestyle=:dash, alpha=0.5)
        end
    end
    
    return fig
end

function CSMakie.pzmap!(ax::Axis, systems::Union{LTISystem, AbstractVector{<:LTISystem}}; 
                        hz=false, kwargs...)
    systems_vec = systems isa AbstractVector ? systems : [systems]
    
    # Add grid at zero
    vlines!(ax, 0, color=:gray, alpha=0.3, linewidth=0.5)
    hlines!(ax, 0, color=:gray, alpha=0.3, linewidth=0.5)
    
    scale_factor = hz ? 1/(2π) : 1
    
    for (i, system) in enumerate(systems_vec)
        p = poles(system) .* scale_factor
        z = tzeros(system) .* scale_factor
        
        # Plot zeros as circles
        if !isempty(z)
            scatter!(ax, real(z), imag(z), 
                    marker=:circle, markersize=15, 
                    color=(:transparent, 0.5), strokewidth=2, 
                    strokecolor=Cycled(i))
        end
        
        # Plot poles as x's
        if !isempty(p)
            scatter!(ax, real(p), imag(p), 
                    marker=:xcross, markersize=15,
                    color=Cycled(i))
        end
        
        # Add unit circle for discrete systems
        if isdiscrete(system) && i == 1  # Only draw once
            θ = range(0, 2π, length=100)
            circle_scale = scale_factor
            lines!(ax, circle_scale .* cos.(θ), circle_scale .* sin.(θ), 
                  color=:gray, linestyle=:dash, alpha=0.5)
        end
    end
    
    return ax
end

# ====== Bodeplot ======
function CSMakie.bodeplot(systems::Union{LTISystem, AbstractVector{<:LTISystem}}, 
                          w=nothing; plotphase=true, unwrap=true, hz=false, 
                          balance=true, adjust_phase_start=true, adaptive=true, kwargs...)
    systems_vec = systems isa AbstractVector ? systems : [systems]
    systems, w = isnothing(w) ? _processfreqplot(Val{:bode}(), systems_vec; adaptive) : 
                                _processfreqplot(Val{:bode}(), systems_vec, w; adaptive)
    
    ws = (hz ? 1/(2π) : 1) .* w
    ny, nu = size(systems[1])
    
    fig = Figure()
    gl = GridLayout(fig[1, 1])
    
    # Create axes grid
    n_rows = plotphase ? 2 : 1
    axes_mag = Matrix{Any}(undef, ny, nu)
    axes_phase = plotphase ? Matrix{Any}(undef, ny, nu) : nothing
    
    for j in 1:nu
        for i in 1:ny
            row_mag = plotphase ? 2i - 1 : i
            
            # Magnitude axis
            ax_mag = Axis(gl[row_mag, j],
                         xscale = log10,
                         yscale = get_yscale_transform(),
                         xlabel = plotphase ? "" : (hz ? "Frequency [Hz]" : "Frequency [rad/s]"),
                         ylabel = i == 1 ? "Magnitude $(ControlSystemsBase._PlotScaleStr)" : "",
                         title = j == 1 && i == 1 ? "Bode Plot" : "")
            axes_mag[i, j] = ax_mag
            
            # Phase axis
            if plotphase
                ax_phase = Axis(gl[2i, j],
                               xscale = log10,
                               xlabel = hz ? "Frequency [Hz]" : "Frequency [rad/s]",
                               ylabel = i == 1 ? "Phase (deg)" : "")
                axes_phase[i, j] = ax_phase
                linkxaxes!(ax_mag, ax_phase)
            end
        end
    end
    
    # Plot data for each system
    for (si, s) in enumerate(systems)
        sbal = balance ? balance_statespace(s)[1] : s
        
        intexcess = 0
        if plotphase && adjust_phase_start && isrational(sbal)
            intexcess = integrator_excess(sbal)
        end
        
        mag, phase = bode(sbal, w; unwrap=false)
        
        if ControlSystemsBase._PlotScale == "dB"
            mag = 20*log10.(mag)
        elseif 0 ∈ mag
            replace!(mag, 0 => NaN)
        end
        
        for j in 1:nu
            for i in 1:ny
                magdata = vec(mag[i, j, :])
                if all(isnan.(magdata) .| isinf.(magdata))
                    continue
                end
                
                # Plot magnitude
                ax_mag = axes_mag[i, j]
                lab = _get_plotlabel(s, i, j)
                
                if adaptive && eltype(magdata) <: AbstractFloat
                    wsi, magi, _ = downsample(ws, magdata, _span(magdata)/500)
                    lines!(ax_mag, wsi, magi, label=lab)
                else
                    lines!(ax_mag, ws, magdata, label=lab)
                end
                
                # Plot phase
                if plotphase
                    phasedata = vec(phase[i, j, :])
                    
                    if adjust_phase_start && isrational(sbal) && intexcess != 0
                        nineties = round(Int, phasedata[1] / 90)
                        phasedata .+= ((90*(-intexcess-nineties)) ÷ 360) * 360
                    end
                    
                    if unwrap
                        phasedata = ControlSystemsBase.unwrap(phasedata .* (π/180)) .* (180/π)
                    end
                    
                    ax_phase = axes_phase[i, j]
                    if adaptive && eltype(phasedata) <: AbstractFloat
                        wsp, phasep, _ = downsample(ws, phasedata, 
                                                   _span(phasedata)/500)
                        lines!(ax_phase, wsp, phasep)
                    else
                        lines!(ax_phase, ws, phasedata)
                    end
                end
            end
        end
    end
    
    return fig
end

# ====== Nyquistplot ======
function CSMakie.nyquistplot(systems::Union{LTISystem, AbstractVector{<:LTISystem}}, 
                             w=nothing; Ms_circles=Float64[], Mt_circles=Float64[], 
                             unit_circle=false, hz=false, critical_point=-1, 
                             balance=true, adaptive=true, kwargs...)
    systems_vec = systems isa AbstractVector ? systems : [systems]
    systems, w = isnothing(w) ? _processfreqplot(Val{:nyquist}(), systems_vec; adaptive) : 
                                _processfreqplot(Val{:nyquist}(), systems_vec, w; adaptive)
    
    ny, nu = size(systems[1])
    
    fig = Figure()
    gl = GridLayout(fig[1, 1])
    
    # Create axes grid
    axes = Matrix{Any}(undef, ny, nu)
    
    for j in 1:nu
        for i in 1:ny
            ax = Axis(gl[i, j],
                     aspect = DataAspect(),
                     xlabel = j == ny ? "Real" : "",
                     ylabel = i == 1 ? "Imaginary" : "",
                     title = i == 1 && j == 1 ? "Nyquist Plot" : "")
            axes[i, j] = ax
            
            # Add grid lines at zero
            vlines!(ax, 0, color=:gray, alpha=0.3, linewidth=0.5)
            hlines!(ax, 0, color=:gray, alpha=0.3, linewidth=0.5)
        end
    end
    
    θ = range(0, 2π, length=100)
    
    for (si, s) in enumerate(systems)
        sbal = balance ? balance_statespace(s)[1] : s
        re_resp, im_resp = nyquist(sbal, w)[1:2]
        
        for j in 1:nu
            for i in 1:ny
                redata = vec(re_resp[i, j, :])
                imdata = vec(im_resp[i, j, :])
                
                ax = axes[i, j]
                lab = _get_plotlabel(s, i, j)
                
                if adaptive
                    indsre = downsample(w, redata, 1/500)[3]
                    indsim = downsample(w, imdata, 1/500)[3]
                    inds = sort!(union(indsre, indsim))
                    lines!(ax, redata[inds], imdata[inds], label=lab)
                else
                    lines!(ax, redata, imdata, label=lab)
                end
                
                # Add circles and critical point for last system
                if si == length(systems)
                    # Ms circles
                    for Ms in Ms_circles
                        lines!(ax, -1 .+ (1/Ms) .* cos.(θ), (1/Ms) .* sin.(θ),
                              color=:gray, linestyle=:dash, alpha=0.5,
                              label="Ms = $(round(Ms, digits=2))")
                    end
                    
                    # Mt circles
                    for Mt in Mt_circles
                        ct = -Mt^2/(Mt^2-1)
                        rt = Mt/(abs(Mt^2-1))
                        lines!(ax, ct .+ rt .* cos.(θ), rt .* sin.(θ),
                              color=:gray, linestyle=:dash, alpha=0.5,
                              label="Mt = $(round(Mt, digits=2))")
                    end
                    
                    # Unit circle
                    if unit_circle
                        lines!(ax, cos.(θ), sin.(θ),
                              color=:gray, linestyle=:dash, alpha=0.5)
                    end
                    
                    # Critical point
                    scatter!(ax, [critical_point], [0],
                            marker=:xcross, markersize=15, color=:red)
                end
            end
        end
    end
    
    return fig
end

# ====== Sigmaplot ======
function CSMakie.sigmaplot(systems::Union{LTISystem, AbstractVector{<:LTISystem}}, 
                           w=nothing; hz=false, balance=true, extrema=false, kwargs...)
    systems_vec = systems isa AbstractVector ? systems : [systems]
    systems, w = isnothing(w) ? _processfreqplot(Val{:sigma}(), systems_vec) : 
                                _processfreqplot(Val{:sigma}(), systems_vec, w)
    
    ws = (hz ? 1/(2π) : 1) .* w
    
    fig = Figure()
    ax = Axis(fig[1,1],
             xscale = log10,
             yscale = get_yscale_transform(),
             title = "Sigma Plot",
             xlabel = hz ? "Frequency [Hz]" : "Frequency [rad/s]",
             ylabel = "Singular Values $(ControlSystemsBase._PlotScaleStr)")
    
    for (si, s) in enumerate(systems)
        sbal = balance ? balance_statespace(s)[1] : s
        sv = sigma(sbal, w)[1]'
        
        if extrema && size(sv, 2) > 2
            sv = sv[:, [1, end]]
        end
        
        if ControlSystemsBase._PlotScale == "dB"
            sv = 20*log10.(sv)
        end
        # Use _to1series to efficiently plot multiple lines
        ws_series, sv_series = _to1series(ws, sv)
        lines!(ax, ws_series, sv_series, color=Cycled(si))
    end
    
    return fig
end

# ====== Marginplot ======
function CSMakie.marginplot(systems::Union{LTISystem, AbstractVector{<:LTISystem}}, 
                            w=nothing; plotphase=true, hz=false, balance=true, 
                            adjust_phase_start=true, adaptive=true, kwargs...)
    systems_vec = systems isa AbstractVector ? systems : [systems]
    systems, w = isnothing(w) ? _processfreqplot(Val{:bode}(), systems_vec; adaptive) : 
                                _processfreqplot(Val{:bode}(), systems_vec, w; adaptive)
    
    ws = (hz ? 1/(2π) : 1) .* w
    ny, nu = size(systems[1])
    
    fig = Figure()
    gl = GridLayout(fig[1, 1])
    
    # Create axes grid
    n_rows = plotphase ? 2 : 1
    axes_mag = Matrix{Any}(undef, ny, nu)
    axes_phase = plotphase ? Matrix{Any}(undef, ny, nu) : nothing
    
    for j in 1:nu
        for i in 1:ny
            row_mag = plotphase ? 2i - 1 : i
            
            # Magnitude axis
            ax_mag = Axis(gl[row_mag, j],
                         xscale = log10,
                         yscale = get_yscale_transform(),
                         xlabel = plotphase ? "" : (hz ? "Frequency [Hz]" : "Frequency [rad/s]"),
                         ylabel = i == 1 ? "Magnitude $_PlotScaleStr[]" : "")
            axes_mag[i, j] = ax_mag
            
            # Phase axis
            if plotphase
                ax_phase = Axis(gl[2i, j],
                               xscale = log10,
                               xlabel = hz ? "Frequency [Hz]" : "Frequency [rad/s]",
                               ylabel = i == 1 ? "Phase (deg)" : "")
                axes_phase[i, j] = ax_phase
                linkxaxes!(ax_mag, ax_phase)
            end
        end
    end
    
    # Plot data for each system
    for (si, s) in enumerate(systems)
        sbal = balance ? balance_statespace(s)[1] : s
        bmag, bphase = bode(sbal, w)
        
        for j in 1:nu
            for i in 1:ny
                wgm, gm, wpm, pm, fullPhase, phasedata = sisomargin(sbal[i,j], w; 
                                                                     full=true, 
                                                                     allMargins=true, 
                                                                     adjust_phase_start)
                
                # Limit number of margins shown
                if length(gm) > 5
                    idx = sortperm(gm)
                    wgm = wgm[idx[1:5]]
                    gm = gm[idx[1:5]]
                end
                if length(pm) > 5
                    idx = sortperm(pm)
                    wpm = wpm[idx[1:5]]
                    pm = pm[idx[1:5]]
                end
                
                # Magnitude plot
                ax_mag = axes_mag[i, j]
                magdata = vec(bmag[i, j, :])
                
                if ControlSystemsBase._PlotScale == "dB"
                    magdata = 20*log10.(magdata)
                    mag_margins = 20 .* log10.(1 ./ gm)
                    oneLine = 0
                else
                    mag_margins = 1 ./ gm
                    oneLine = 1
                end
                
                # Plot magnitude response
                if adaptive
                    wsi, magi, _ = downsample(ws, magdata, 
                                            _span(magdata)/500)
                    lines!(ax_mag, wsi, magi, label=_get_plotlabel(s, i, j))
                else
                    lines!(ax_mag, ws, magdata, label=_get_plotlabel(s, i, j))
                end
                
                # Unity gain line
                hlines!(ax_mag, oneLine, color=:gray, linestyle=:dash, alpha=0.5)
                
                # Gain margins - draw vertical lines from unity to margin value
                wgm_display = hz ? wgm ./ (2π) : wgm
                for k in eachindex(wgm_display)
                    # Draw vertical line from unity gain to the margin value
                    lines!(ax_mag, [wgm_display[k], wgm_display[k]], 
                          [oneLine, mag_margins[k]], color=:red, linewidth=2)
                end
                
                # Add title with margins info
                if si == length(systems)
                    gm_str = "Gm: [" * join([Printf.@sprintf("%.2g", g) for g in gm], ", ") * "]"
                    wgm_str = "ωgm: [" * join([Printf.@sprintf("%.2g", w) for w in wgm_display], ", ") * "]"
                    ax_mag.title = gm_str * " " * wgm_str
                end
                
                # Phase plot
                if plotphase
                    ax_phase = axes_phase[i, j]
                    
                    if adaptive
                        wsp, phasep, _ = downsample(ws, phasedata, 
                                                  _span(phasedata)/500)
                        lines!(ax_phase, wsp, phasep)
                    else
                        lines!(ax_phase, ws, phasedata)
                    end
                    
                    # Phase margin lines
                    wpm_display = hz ? wpm ./ (2π) : wpm
                    
                    # Draw horizontal lines at phase margin crossings
                    for k in 1:length(pm)
                        phase_line = fullPhase[k] - pm[k]
                        hlines!(ax_phase, phase_line, color=:gray, linestyle=:dash, alpha=0.5)
                    end
                    
                    # Draw vertical lines showing the phase margins
                    for k in 1:length(pm)
                        lines!(ax_phase, [wpm_display[k], wpm_display[k]], 
                              [fullPhase[k], fullPhase[k] - pm[k]], color=:blue, linewidth=2)
                    end
                    
                    # Add title with phase margins info
                    if si == length(systems)
                        pm_str = "Pm: [" * join([Printf.@sprintf("%.2g°", p) for p in pm], ", ") * "]"
                        wpm_str = "ωpm: [" * join([Printf.@sprintf("%.2g", w) for w in wpm_display], ", ") * "]"
                        ax_phase.title = pm_str * " " * wpm_str
                    end
                end
            end
        end
    end
    
    return fig
end

# ====== Root Locus Plot ======
function CSMakie.rlocusplot(P::LTISystem, K=500; output=false, kwargs...)
    # Compute root locus
    result = rlocus(P, K; output=output)
    roots, Z, K_vals = result.roots, result.Z, result.K
    
    array_K = eltype(K_vals) <: AbstractArray
    
    fig = Figure()
    ax = Axis(fig[1,1],
             aspect = DataAspect(),
             title = "Root Locus",
             xlabel = "Re(roots)",
             ylabel = "Im(roots)")
    
    # Add grid at zero
    vlines!(ax, 0, color=:gray, alpha=0.3, linewidth=0.5)
    hlines!(ax, 0, color=:gray, alpha=0.3, linewidth=0.5)
    
    redata = real.(roots)
    imdata = imag.(roots)
    
    # Plot root locus branches
    for i in 1:size(redata, 2)
        lines!(ax, redata[:, i], imdata[:, i], color=:blue)
    end
    
    # Plot zeros
    if !isempty(Z)
        scatter!(ax, real.(Z), imag.(Z), 
                marker=:circle, markersize=10, 
                color=(:transparent, 0.5), strokewidth=2, 
                strokecolor=:green, label="Zeros")
    end
    
    # Plot open-loop poles
    scatter!(ax, redata[1, :], imdata[1, :], 
            marker=:xcross, markersize=10, 
            color=:red, label="Open-loop poles")
    
    # Plot closed-loop poles for matrix K
    if array_K
        scatter!(ax, redata[end, :], imdata[end, :], 
                marker=:diamond, markersize=10, 
                color=:purple, label="Closed-loop poles")
    end
    
    return fig
end

# ====== RGA Plot ======
function CSMakie.rgaplot(systems::Union{LTISystem, AbstractVector{<:LTISystem}}, 
                         w=nothing; hz=false, balance=true, kwargs...)
    systems_vec = systems isa AbstractVector ? systems : [systems]
    systems, w = isnothing(w) ? _processfreqplot(Val{:sigma}(), systems_vec) : 
                                _processfreqplot(Val{:sigma}(), systems_vec, w)
    
    ws = (hz ? 1/(2π) : 1) .* w
    
    fig = Figure()
    ax = Axis(fig[1,1],
             xscale = log10,
             title = "RGA Plot",
             xlabel = hz ? "Frequency [Hz]" : "Frequency [rad/s]",
             ylabel = "Element magnitudes")
    
    for (si, s) in enumerate(systems)
        sbal = balance ? balance_statespace(s)[1] : s
        rga = abs.(relative_gain_array(sbal, w))
        
        for j in 1:size(rga, 1)
            for i in 1:size(rga, 2)
                label = "System $si, from $i to $j"
                lines!(ax, ws, rga[j, i, :], label=label)
            end
        end
    end
    
    Legend(fig[1, 2], ax)
    
    return fig
end

# ====== LeadLinkCurve ======
function CSMakie.leadlinkcurve(start=1; kwargs...)
    N = range(start, stop=10, length=50)
    dph = map(Ni -> 180/π*atan(sqrt(Ni)) - atan(1/sqrt(Ni)), N)
    
    fig = Figure()
    ax = Axis(fig[1,1],
             xlabel = "N",
             ylabel = "Phase advance [deg]",
             title = "Lead Link Phase Advance")
    
    lines!(ax, N, dph, color=:blue, linewidth=2)
    
    return fig
end

# ====== Nicholsplot ======
# Note: This is a simplified version. The full implementation would require
# porting all the complex grid calculations from the Plots version
function CSMakie.nicholsplot(systems::Union{LTISystem, AbstractVector{<:LTISystem}}, 
                             w=nothing; text=true, 
                             Gains=[12, 6, 3, 1, 0.5, -0.5, -1, -3, -6, -10, -20, -40, -60],
                             pInc=30, kwargs...)
    systems_vec = systems isa AbstractVector ? systems : [systems]
    systems, w = isnothing(w) ? _processfreqplot(Val{:nyquist}(), systems_vec) : 
                                _processfreqplot(Val{:nyquist}(), systems_vec, w)
    
    fig = Figure()
    ax = Axis(fig[1,1],
             title = "Nichols Chart",
             xlabel = "Phase [deg]",
             ylabel = "Magnitude [dB]")
    
    # Plot systems
    for (sysi, s) in enumerate(systems)
        ℜresp, ℑresp = nyquist(s, w)[1:2]
        ℜdata = dropdims(ℜresp, dims=(1,2))
        ℑdata = dropdims(ℑresp, dims=(1,2))
        mag = 20*log10.(sqrt.(ℜdata.^2 + ℑdata.^2))
        angles = 180/π*angle.(ℜdata .+ im*ℑdata)
        
        lines!(ax, angles, mag, linewidth=2, label="System $sysi")
    end
    
    # Note: Full Nichols grid implementation would go here
    # This is a placeholder for the complex grid calculations
    
    return fig
end

# ====== Recipes for Types (SimResult, StepInfo) ======
# These use Makie's standard recipe system for plotting types

@recipe(SimResultPlot, r) do scene
    Attributes(
        plotu = false,
        plotx = false,
        ploty = true,
    )
end

function Makie.plot!(srp::SimResultPlot)
    r = srp.r[]
    
    # Handle vector of SimResults
    if r isa AbstractVector
        for res in r
            plot!(current_axis(), res; 
                  plotu=srp.plotu[], plotx=srp.plotx[], ploty=srp.ploty[])
        end
        return srp
    end
    
    plotu = srp.plotu[]
    plotx = srp.plotx[]
    ploty = srp.ploty[]
    
    ny, nu, nx = r.ny, r.nu, r.nx
    t = r.t
    n_series = size(r.y, 3)
    
    # Calculate number of subplots needed
    nplots = ploty ? ny : 0
    plotu && (nplots += nu)
    plotx && (nplots += nx)
    
    # Get current figure or create new one
    fig = current_figure()
    
    # Create grid layout for subplots
    gl = GridLayout(fig[1, 1])
    
    plotind = 1
    axes = []
    
    # Plot outputs
    if ploty
        for i in 1:ny
            ax = Axis(gl[plotind, 1], 
                     xlabel = plotind == nplots ? "Time (s)" : "",
                     ylabel = output_names(r.sys)[i])
            push!(axes, ax)
            
            for ms in 1:n_series
                label = n_series > 1 ? "From $(input_names(r.sys)[ms])" : ""
                if iscontinuous(r.sys)
                    lines!(ax, t, r.y[i, :, ms], label=label)
                else
                    stairs!(ax, t, r.y[i, :, ms], label=label, step=:post)
                end
            end
            plotind += 1
        end
    end
    
    # Plot inputs
    if plotu
        for i in 1:nu
            ax = Axis(gl[plotind, 1],
                     xlabel = plotind == nplots ? "Time (s)" : "",
                     ylabel = input_names(r.sys)[i])
            push!(axes, ax)
            
            if iscontinuous(r.sys)
                lines!(ax, t, r.u[i, :])
            else
                stairs!(ax, t, r.u[i, :], step=:post)
            end
            plotind += 1
        end
    end
    
    # Plot states
    if plotx
        for i in 1:nx
            ax = Axis(gl[plotind, 1],
                     xlabel = plotind == nplots ? "Time (s)" : "",
                     ylabel = state_names(r.sys)[i])
            push!(axes, ax)
            
            if n_series > 1
                for ms in 1:n_series
                    if iscontinuous(r.sys)
                        lines!(ax, t, r.x[i, :, ms])
                    else
                        stairs!(ax, t, r.x[i, :, ms], step=:post)
                    end
                end
            else
                if iscontinuous(r.sys)
                    lines!(ax, t, r.x[i, :])
                else
                    stairs!(ax, t, r.x[i, :], step=:post)
                end
            end
            plotind += 1
        end
    end
    
    # Link x-axes
    if length(axes) > 1
        linkxaxes!(axes...)
    end
    
    srp
end

@recipe(StepInfoPlot, si) do scene
    Attributes()
end

function Makie.plot!(sip::StepInfoPlot)
    si = sip.si[]
    
    ax = current_axis()
    
    # Plot the step response
    t = si.res.t
    y = si.res.y[1, :]
    
    lines!(ax, t, y, color=:blue, linewidth=2)
    
    # Final value
    hlines!(ax, si.yf, color=:black, linestyle=:dash, 
           label=Printf.@sprintf("Final value: %.3f", si.yf))
    
    # Rise time
    lines!(ax, t[si.i10:si.i90], y[si.i10:si.i90], 
          color=:green, linewidth=3,
          label=Printf.@sprintf("Rise time: %.3f", si.risetime))
    
    vlines!(ax, [t[si.i10], t[si.i90]], color=:green, linestyle=:dash, alpha=0.5)
    
    # Peak and overshoot
    scatter!(ax, [si.peaktime], [si.peak], color=:red, markersize=10,
            label=Printf.@sprintf("Peak: %.3f, Overshoot: %.1f%%", si.peak, si.overshoot))
    lines!(ax, [si.peaktime, si.peaktime], [si.y0, si.peak], color=:red, alpha=0.5)
    
    # Settling time
    scatter!(ax, [si.settlingtime], [y[si.settlingtimeind]], color=:purple, markersize=10,
            label=Printf.@sprintf("Settling time: %.3f", si.settlingtime))
    lines!(ax, [si.settlingtime, si.settlingtime], [si.y0, y[si.settlingtimeind]], 
          color=:purple, alpha=0.5)
    
    # Settling threshold bands
    hlines!(ax, [si.yf - si.stepsize*si.settling_th, si.yf + si.stepsize*si.settling_th],
           color=:purple, linestyle=:dash, alpha=0.3,
           label=Printf.@sprintf("Settling threshold: %.1f%%", 100*si.settling_th))
    
    # Undershoot if present
    if si.undershoot != 0
        t_under = t[si.lowerpeakind]
        scatter!(ax, [t_under], [si.lowerpeak], color=:orange, markersize=10,
                label=Printf.@sprintf("Undershoot: %.1f%%", si.undershoot))
        lines!(ax, [t_under, t_under], [si.y0, si.lowerpeak], color=:orange, alpha=0.5)
    end
    
    ax.xlabel = "Time (s)"
    ax.ylabel = "Output"
    ax.title = "Step Response Analysis"
    
    axislegend(ax, position=:rt)
    
    sip
end

end # module