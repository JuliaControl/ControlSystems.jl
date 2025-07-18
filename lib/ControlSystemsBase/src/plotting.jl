export bodeplot, nyquistplot, sigmaplot, marginplot, setPlotScale, gangoffour, gangoffourplot, gangofseven, pzmap, pzmap!, nicholsplot

_PlotScale = "log10"
_PlotScaleFunc = :log10
_PlotScaleStr = ""



"""
    setPlotScale(str)

Set the default scale of magnitude in `bodeplot` and `sigmaplot`.
`str` should be either `"dB"` or `"log10"`. The default scale if none is chosen is `"log10"`.
"""
function setPlotScale(str::AbstractString)
    if str == "dB"
        plotSettings = (str, :identity, "(dB)")
    elseif str == "log10"
        plotSettings = (str, :log10, "")
    else
        error("Scale must be set to either \"dB\" or \"log10\"")
    end
    global _PlotScale, _PlotScaleFunc, _PlotScaleStr
    _PlotScale, _PlotScaleFunc, _PlotScaleStr = plotSettings
end

# """
# Get attributes from xlims or ylims
# default to extrema(wmag) if xlims/ylims not defined or empty
# """
# function getlims(xylims, plotattributes, wmag)
#     lims = get(plotattributes, xylims, extrema(wmag))
#     if !isa(lims, Tuple{<:Number, <:Number}) # If x/ylims not supplied as empty
#         lims = extrema(wmag)
#     end
#     if !isempty(get_serieslist(plotattributes))
#         subplot = get(plotattributes, :subplot, 0)
#         (subplot == 0 || (subplot isa Array)) && (return lims)
#         se = seriesextrema(xylims, plotattributes, subplot)
#         lims = extremareducer(lims, se)
#     end
#     lims
# end

# get_serieslist(plotattributes) = plotattributes[:plot_object].series_list
# get_serieslist(plotattributes, subplot) = plotattributes[:plot_object].subplots[subplot].series_list

# function seriesextrema(xylims, plotattributes, subplot)
#     serieslist = get_serieslist(plotattributes, subplot)
#     isempty(serieslist) && (return (Inf, -Inf))
#     sym = xylims === :xlims ? :x : :y
#     mapreduce(extremareducer, serieslist) do series
#         extrema(series[sym])
#     end
# end
# extremareducer(x,y) = (min(x[1],y[1]),max(x[2],y[2]))

# function getPhaseTicks(x, minmax)
#     minx, maxx =  minmax
#     min               = ceil(minx/90)
#     max               = floor(maxx/90)
#     if max-min < 5
#         ## If we span less than a full rotation 45° steps are ok
#         major = ((min-0.5):0.5:(max+0.5)).*90
#     else
#         ## Create additional 45° before/behind first/last plot
#         ## this helps identifying at the edges.
#         major = [(min-0.5);min:max;(max+0.5)].*90
#     end
#     majorText = ["$(round(Int64,i))" for i = major]

#     return major, majorText

# end

# function getLogTicks(x, minmax)
#     minx, maxx =  minmax
#     major_minor_limit = 6
#     minor_text_limit  = 8
#     min               = minx <= 0 ? minimum(x) : ceil(log10(minx))
#     max               = floor(log10(maxx))
#     major             = exp10.(min:max)

#     majorText = ["10^{$(round(Int64,i))}" for i = min:max]

#     if max - min < major_minor_limit
#         minor     = [j*exp10(i) for i = (min-1):(max+1) for j = 2:9]
#         minorText = ["$j*10^{$(round(Int64,i))}" for i = (min-1):(max+1) for j = 2:9]
#         ind       = findall(minx .<= minor .<= maxx)
#         minor     = minor[ind]
#         minorText = minorText[ind]
#         if length(minor) > minor_text_limit
#             minorText = [" " for t in minorText]#fill!(minorText, " ")
#         end
#         perm = sortperm([major; minor])
#         return [major; minor][perm], [majorText; minorText][perm]

#     else
#         return major, majorText
#     end
# end


# This will be called on plot(lsim(sys, args...))
@recipe function simresultplot(r::SimResult; plotu=false, plotx=false, ploty=true, input_names=ControlSystemsBase.input_names(r.sys), output_names=ControlSystemsBase.output_names(r.sys), state_names=ControlSystemsBase.state_names(r.sys))
    ny, nu, nx, sys = r.ny, r.nu, r.nx, r.sys
    t = r.t
    n_series = size(r.y, 3) # step and impulse produce multiple results
    nplots = ploty ? ny : 0
    plotu && (nplots += nu)
    plotx && (nplots += nx)
    layout --> (nplots, 1)
    seriestype --> (iscontinuous(r.sys) ? :path : :steppost)
    plotind = 1
    if ploty
        for ms in 1:n_series
            for i=1:ny
                ytext = output_names[i]
                @series begin
                    xguide  --> "Time (s)"
                    yguide  --> ytext
                    label   --> (n_series > 1 ? "From $(input_names[ms])" : "")
                    subplot --> i
                    t,  r.y[i, :, ms]
                end
            end
        end 
        plotind = ny+1
    end 
    if plotu # bug in recipe system, can't use `plotu || return`
        for i=1:nu
            utext = input_names[i]
            @series begin
                xguide  --> "Time (s)"
                yguide  --> utext
                subplot --> plotind
                label --> ""
                t,  r.u[i, :]
            end
            plotind += 1
        end
    end
    if plotx
        for i=1:nx
            xtext = state_names[i]
            @series begin
                xguide  --> "Time (s)"
                yguide  --> xtext
                subplot --> plotind
                label --> ""
                t,  (n_series > 1 ? r.x[i, :, :] : r.x[i, :])
            end
            plotind += 1
        end
    end
end

@recipe function simresultplot(r::AbstractVector{<:SimResult})
    for r in r
        @series begin
            r
        end
    end
end

@recipe function stepinfoplot(si::StepInfo)
    @series begin
        color --> 1
        si.res
    end
    linestyle --> :dash
    @series begin
        color --> 1
        seriestype := :hline
        label := @sprintf("Final value: %.3f", si.yf)
        [si.yf]
    end
    @series begin
        linestyle := :solid
        linewidth --> 2
        color --> 2
        label := @sprintf("Rise time: %.3f", si.risetime)
        si.res.t[si.i10:si.i90], si.res.y[1, si.i10:si.i90]
    end
    @series begin
        color --> 2
        seriestype := :vline
        label := @sprintf("Rise time threshold: %.1f%%-%.1f%%", 100si.risetime_th[1], 100si.risetime_th[2])
        [si.res.t[si.i10], si.res.t[si.i90]]
    end
    @series begin
        color --> 3
        label := @sprintf("Peak: %.3f Overshoot: %.1f%%", si.peak, si.overshoot)
        markershape --> [:none, :circle]
        si.peaktime*ones(2), [si.y0, si.peak]
    end
    @series begin
        color --> 4
        markershape --> [:none, :circle]
        label := @sprintf("Settling time: %.3f", si.settlingtime)
        si.settlingtime*ones(2), [si.y0, si.res.y[1, si.settlingtimeind]]
    end
    @series begin
        color --> 4
        seriestype := :hline
        label := @sprintf("Settling threshold: %.1f%%", 100si.settling_th)
        [si.yf-si.stepsize*si.settling_th, si.yf+si.stepsize*si.settling_th]
    end
    if si.undershoot != 0
        @series begin
            color --> 5
            label := @sprintf("Undershoot: %.1f%%", si.undershoot)
            markershape --> [:none, :circle]
            t = si.res.t[si.lowerpeakind]
            t*ones(2), [si.y0, si.lowerpeak]
        end
    end

end

"""
    _processfreqplot(plottype, system::LTISystem, [w])
    _processfreqplot(plottype, system::AbstractVector{<:LTISystem}, [w])

Calculate default frequency vector and put system in array of not already array.
`plottype` is one of `Val{:bode}, Val{:nyquist}, ...`
for which `_default_freq_vector` is defined.
Check that system dimensions are compatible.
"""
_processfreqplot(plottype, system::LTISystem, args...; kwargs...) =
    _processfreqplot(plottype, [system], args...; kwargs...)
# Catch when system is not vector, with and without frequency input

# Catch correct form
_processfreqplot(plottype, systems::AbstractVector{<:LTISystem}; adaptive=false) =
    _processfreqplot(plottype, systems, _default_freq_vector(systems, plottype; adaptive))

function _processfreqplot(plottype, systems::AbstractVector{<:LTISystem}, w; kwargs...)

    if !_same_io_dims(systems...)
        error("All systems must have the same input/output dimensions")
    end
    return systems, w
end


@userplot Bodeplot
## FREQUENCY PLOTS ##
"""
    fig = bodeplot(sys, args...)
    bodeplot(LTISystem[sys1, sys2...], args...; plotphase=true, balance = true, kwargs...)

Create a Bode plot of the `LTISystem`(s). A frequency vector `w` can be
optionally provided. To change the Magnitude scale see [`setPlotScale`](@ref). The default magnitude scale is "log10" (absolute scale).
                                            
- If `hz=true`, the plot x-axis will be displayed in Hertz, the input frequency vector is still treated as rad/s.
- `balance`: Call [`balance_statespace`](@ref) on the system before plotting.
- `adjust_phase_start`: If true, the phase will be adjusted so that it starts at -90*intexcess degrees, where `intexcess` is the integrator excess of the system.
- `adaptive`: If true, an adaptive frequency grid is used in order to keep the number of plotted points low, while resolving features in the frequency response well. If a manually provided frequency vector is used, this may be downsampled before plotting.

`kwargs` is sent as argument to RecipesBase.plot.
"""
bodeplot

function _get_plotlabel(s, i, j)
    sys_name = system_name(s)
    if !isempty(sys_name)
        sys_name = sys_name * ": "
    end
    u_names = input_names(s)
    y_names = output_names(s)
    default_names = all(match(r"^u\(?\d*\)?$", name) !== nothing for name in u_names) &&
        all(match(r"^y\(?\d*\)?$", name) !== nothing for name in y_names)
    if default_names && isempty(sys_name)
        nothing
    else # It's enough that the system has either a system name or signal names
        "$(sys_name)$(u_names[j]) → $(y_names[i])"
    end
end

_span(vec) = -(reverse(extrema(vec))...)

@recipe function bodeplot(p::Bodeplot; plotphase=true, ylimsphase=(), unwrap=true, hz=false, balance=true, adjust_phase_start=true, adaptive=true)
    systems, w = _processfreqplot(Val{:bode}(), p.args...; adaptive)
    ws = (hz ? 1/(2π) : 1) .* w
    ny, nu = size(systems[1])
    s2i(i,j) = LinearIndices((nu,(plotphase ? 2 : 1)*ny))[j,i]
    layout --> ((plotphase ? 2 : 1)*ny, nu)
    nw = length(w)
    # xticks --> getLogTicks(ws, getlims(:xlims, plotattributes, ws))
    grid   --> true

    for (si,s) = enumerate(systems)
        if balance
            sbal = balance_statespace(s)[1]
        else
            sbal = s
        end
        if plotphase && adjust_phase_start && isrational(sbal)
            intexcess = integrator_excess(sbal)
        end
        mag, phase = bode(sbal, w; unwrap=false)
        if _PlotScale == "dB" # Set by setPlotScale(str) globally
            mag = 20*log10.(mag)
        elseif 0 ∈ mag
            replace!(mag, 0 => -Inf) # To prevent plot crashing when some magnitude is exactly zero
        end

        xlab = plotphase ? "" : (hz ? "Frequency [Hz]" : "Frequency [rad/s]")
        group_ind = 0
        for j=1:nu
            for i=1:ny
                group_ind += 1
                magdata = vec(mag[i, j, :])
                if all(magdata .== -Inf)
                    # 0 system, don't plot anything
                    continue
                end
                phasedata = vec(phase[i, j, :])
                local inds
                @series begin
                    yscale    --> _PlotScaleFunc
                    xscale    --> :log10
                    # if _PlotScale != "dB"
                    #     yticks    --> getLogTicks(magdata, getlims(:ylims, plotattributes, magdata))
                    # end
                    xguide    --> xlab
                    yguide    --> "Magnitude $_PlotScaleStr"
                    subplot   --> min(s2i((plotphase ? (2i-1) : i),j), prod(plotattributes[:layout]))
                    lab = _get_plotlabel(s, i, j)
                    if lab !== nothing
                        label --> lab
                    end
                    group     --> group_ind
                    if adaptive
                        lmag = _PlotScale == "dB" ? magdata : log.(magdata)
                        wsi, _, inds = downsample(ws, lmag, _span(lmag)/500)
                        wsi, magdata[inds]
                    else
                        ws, magdata
                    end
                end
                plotphase || continue

                if adjust_phase_start == true && isrational(sbal)
                    if intexcess != 0
                        # Snap phase so that it starts at -90*intexcess
                        nineties = round(Int, phasedata[1] / 90)
                        phasedata .+= ((90*(-intexcess-nineties)) ÷ 360) * 360
                    end
                end

                if eltype(phasedata) <: AbstractFloat
                    link --> :x # To guard agains https://github.com/JuliaPlots/Plots.jl/issues/5092 when using uncertain number systems
                end

                @series begin
                    xscale    --> :log10
                    # ylims      := ylimsphase
                    # yticks    --> yphaseticks
                    yguide    --> "Phase (deg)"
                    subplot   --> s2i(2i,j)
                    xguide    --> (hz ? "Frequency [Hz]" : "Frequency [rad/s]")
                    label     --> ""
                    group     --> group_ind
                    phasedata = unwrap ? ControlSystemsBase.unwrap(phasedata.*(pi/180)).*(180/pi) : phasedata
                    if adaptive && eltype(phasedata) <: AbstractFloat # To guard agains https://github.com/JuliaPlots/Plots.jl/issues/5092 when using uncertain number systems
                        downsample(ws, phasedata, _span(phasedata)/500)[1:2]
                    elseif adaptive
                        ws[inds], phasedata[inds]
                    else
                        ws, phasedata
                    end

                end
            end
        end
    end
end

@recipe function f(::Type{Val{:bodemag}}, x, y, z)
    w = x
    magdata = y
    seriestype := :path
    primary --> false
    grid   --> true
    yscale --> _PlotScaleFunc
    xscale --> :log10
    yguide --> "Magnitude $_PlotScaleStr"
    x := w
    y := magdata
    ()
end
@recipe function f(::Type{Val{:bodephase}}, x, y, z)
    w = x
    phasedata = y
    seriestype := :path
    primary --> false
    grid   --> true
    xscale --> :log10
    yguide --> "Phase (deg)"
    xguide --> "Frequency (rad/s)"
    x := w
    y := phasedata
    ()
end



@userplot Nyquistplot
"""
    fig = nyquistplot(sys;                Ms_circles=Float64[], Mt_circles=Float64[], unit_circle=false, hz=false, critical_point=-1, kwargs...)
    nyquistplot(LTISystem[sys1, sys2...]; Ms_circles=Float64[], Mt_circles=Float64[], unit_circle=false, hz=false, critical_point=-1, kwargs...)

Create a Nyquist plot of the `LTISystem`(s). A frequency vector `w` can be
optionally provided.

- `unit_circle`: if the unit circle should be displayed. The Nyquist curve crosses the unit circle at the gain crossover frequency.
- `Ms_circles`: draw circles corresponding to given levels of sensitivity (circles around -1 with  radii `1/Ms`). `Ms_circles` can be supplied as a number or a vector of numbers. A design staying outside such a circle has a phase margin of at least `2asin(1/(2Ms))` rad and a gain margin of at least `Ms/(Ms-1)`. See also [`margin_bounds`](@ref), [`Ms_from_phase_margin`](@ref) and [`Ms_from_gain_margin`](@ref).
- `Mt_circles`: draw circles corresponding to given levels of complementary sensitivity. `Mt_circles` can be supplied as a number or a vector of numbers.
- `critical_point`: point on real axis to mark as critical for encirclements
- If `hz=true`, the hover information will be displayed in Hertz, the input frequency vector is still treated as rad/s.
- `balance`: Call [`balance_statespace`](@ref) on the system before plotting.

`kwargs` is sent as argument to plot.
"""
nyquistplot
@recipe function nyquistplot(p::Nyquistplot; Ms_circles=Float64[], Mt_circles=Float64[], unit_circle=false, hz=false, critical_point=-1, balance=true, adaptive=true)
    systems, w = _processfreqplot(Val{:nyquist}(), p.args...; adaptive)
    ny, nu = size(systems[1])
    nw = length(w)
    layout --> (ny,nu)
    framestyle --> :zerolines
    s2i(i,j) = LinearIndices((nu,ny))[j,i]
    θ = range(0, stop=2π, length=100)
    S, C = sin.(θ), cos.(θ)
    for (si,s) = enumerate(systems)
        if balance
            s = balance_statespace(s)[1]
        end
        re_resp, im_resp = nyquist(s, w)[1:2]
        for j=1:nu
            for i=1:ny
                redata = re_resp[i, j, :]
                imdata = im_resp[i, j, :]
                mask = @. (-20 ≤ imdata ≤ 20) & (-20 ≤ redata ≤ 20)
                ylims --> (min(minimum(imdata[mask]),-1.05), max(maximum(imdata[mask]),1.05))
                xlims --> (min(minimum(redata[mask]),-1.05), max(maximum(redata[mask]),1.05))
                @series begin
                    subplot --> s2i(i,j)
                    lab = _get_plotlabel(s, i, j)
                    if lab !== nothing
                        label --> lab
                    end
                    hover_data = [hz ? Printf.@sprintf("f = %.3g", w/2π) : Printf.@sprintf("ω = %.3g", w) for w in w]
                    if adaptive
                        indsre = downsample(w, redata, 1/500)[3]
                        indsim = downsample(w, imdata, 1/500)[3]
                        inds = sort!(union(indsre, indsim))
                        hover --> hover_data[inds]
                        redata[inds], imdata[inds]
                    else
                        hover --> hover_data
                        redata, imdata
                    end
                end                
                
                if si == length(systems)
                    for Ms in Ms_circles
                        @series begin
                            subplot --> s2i(i,j)
                            primary := false
                            linestyle := :dash
                            linecolor := :gray
                            seriestype := :path
                            markershape := :none
                            label := "Ms = $(round(Ms, digits=2))"
                            (-1 .+ (1/Ms) * C, (1/Ms) * S)
                        end
                    end 
                    for Mt in Mt_circles
                        @series begin
                            subplot --> s2i(i,j)
                            primary := false
                            linestyle := :dash
                            linecolor := :gray
                            seriestype := :path
                            markershape := :none
                            label := "Mt = $(round(Mt, digits=2))"
                            ct = -Mt^2/(Mt^2-1) # Mt center
                            rt = Mt/(Mt^2-1)    # Mt radius
                            ct.+rt.*C, rt.*S
                        end
                    end                
                    if unit_circle 
                        @series begin
                            subplot --> s2i(i,j)
                            primary := false
                            linestyle := :dash
                            linecolor := :gray
                            seriestype := :path
                            markershape := :none
                            (C, S)
                        end
                    end
                    @series begin # Mark the critical point
                        # Title and yguide must be here in the last series for the result to be correct
                        title --> "Nyquist plot from: $(input_names(s, j))"
                        yguide --> "To: $(output_names(s, i))"
                        subplot --> s2i(i,j)
                        primary := false
                        markershape := :xcross
                        seriescolor := :red
                        markersize := 5
                        seriestype := :scatter
                        [critical_point], [0]
                    end
                end
                
            end
        end
    end
end


@userplot Nicholsplot

"""
    fig = nicholsplot{T<:LTISystem}(systems::Vector{T}, w::AbstractVector; kwargs...)

Create a Nichols plot of the `LTISystem`(s). A frequency vector `w` can be
optionally provided.

Keyword arguments:
```
text = true
Gains = [12, 6, 3, 1, 0.5, -0.5, -1, -3, -6, -10, -20, -40, -60]
pInc = 30
sat = 0.4
val = 0.85
fontsize = 10
```

`pInc` determines the increment in degrees between phase lines.

`sat` ∈ [0,1] determines the saturation of the gain lines

`val` ∈ [0,1] determines the brightness of the gain lines

Additional keyword arguments are sent to the function plotting the systems and can be
used to specify colors, line styles etc. using regular RecipesBase.jl syntax

This function is based on code subject to the two-clause BSD licence
Copyright 2011 Will Robertson
Copyright 2011 Philipp Allgeuer

"""
nicholsplot
@recipe function nicholsplot(p::Nicholsplot;
    text     = true,
    Gains    = [12, 6, 3, 1, 0.5, -0.5, -1, -3, -6, -10, -20, -40, -60],
    pInc     = 30,
    sat      = 0.4,
    val      = 0.85,
    fontsize = 10)

    plots_id = Base.PkgId(UUID("91a5bcdd-55d7-5caf-9e0b-520d859cae80"), "Plots")
    haskey(Base.loaded_modules, plots_id) || error("Call using Plots before calling this function")
    Plots = Base.loaded_modules[plots_id]
    Colors = Plots.Colors

    systems, w = _processfreqplot(Val{:nyquist}(), p.args...)
    ny, nu = size(systems[1])

    if isdiscrete(systems[1])
        w_nyquist = 2π/systems[1].Ts
        w = w[w.<= w_nyquist]
    end
    nw = length(w)

    # Gain circle functions
    angle(x)        = unwrap(atan.(imag.(x),real.(x)))
    RadM(m)         = @. abs(m/(m^2-1))
    CentreM(m)      = @. m^2/(1-m^2)
    Ny(mdb,t)       = @. CentreM(10^(mdb/20))+RadM(10^(mdb/20)).*(cosd(t)+im.*sind(t))
    Niϕ(mdb,t)      = @. rad2deg((angle(Ny(mdb,t))))
    Ni_Ga(mdb,t)    = @. 20 * log10(abs(Ny(mdb,t)))

    # Phase circle functions
    Radϕ(ϕ)         = @. 1 / (2 * abs(sind(ϕ)))
    Nyℜ(ϕ,t)        = @. -0.5+Radϕ(ϕ)*cosd(t+mod(ϕ,180)-90)
    Nyℑ(ϕ,t)        = @. 1 / (2 .* tand(ϕ))+Radϕ(ϕ)*sind(t+mod(ϕ,180)-90)
    Niϕϕ(ϕ,t)       = @. rad2deg((angle(Nyℜ(ϕ,t)+im*Nyℑ(ϕ,t))))+360*(round(ϕ/360,RoundToZero)+(0t<0))
    Ni_Gaϕ(ϕ,t)     = @. 20 * log10(abs(Nyℜ(ϕ,t)+im*Nyℑ(ϕ,t)))
    Ni_La(ϕ)        = @. 0.090*10^(ϕ/60)
    getColor(mdb)   = convert(Colors.RGB,Colors.HSV(360*((mdb-minimum(Gains))/(maximum(Gains)-minimum(Gains)))^1.5,sat,val))

    megaangles      = vcat(map(s -> 180/π*angle(vec(freqresp(s, w))), systems)...)
    filter!(x-> !isnan(x), megaangles)
    extremeangles = extrema(megaangles)
    extremeangles = floor(extremeangles[1]/180)*180, ceil(extremeangles[2]/180)*180
    PCyc            = Set{Int}(floor.(Int,megaangles/360)) |> collect |> sort
    # PCyc            = extremeangles[1]:pInc:extremeangles[2]

    # yticks := (Float64[],String[])
    yguide --> "Open-loop gain [dB]"
    xguide --> "Open-loop phase [deg]"

    #  Gain circles
    for k=Gains
        ϕVals   =Niϕ(k,-180:1:180)
        GVals   =Ni_Ga(k,-180:1:180)
        for l in PCyc
            @series begin
                linewidth := 1
                linecolor := getColor(k)
                grid --> false
                if text
                    offset  = (l+1)
                    TextX   = Niϕ(k,210) .+offset
                    TextY   = Ni_Ga(k,210)
                    annotations := (TextX,TextY,("$(string(k)) dB"))
                end
                ϕVals .+ 360(l+1),GVals
            end
        end
    end

    #  Phase circles
    PCycbottom = (PCyc[1] < 0 ? PCyc[1]*360 : (PCyc[1]-1)*360)
    PCyctop = (PCyc[end] < 0 ? (PCyc[end]+1)*360 : (PCyc[end])*360)

    Phi=(PCycbottom):pInc:(PCyctop)
    T1 = 10.0 .^range(-4,stop=log10(180), length=300)
    T2 = [T1; 360 .- reverse(T1,dims=1)]

    for k=(Phi .+ 180)
        if abs(sind(k))<1e-3
            @series begin
                linewidth := 1
                linecolor := Colors.RGB(0.75*[1, 1, 1]...)
                [k,k],[-110,25]
            end
            if cosd(5)>0
                TextX=k
                TextY=1
            else
                TextX=k
                TextY=-46.5
            end
        else
            @series begin
                linewidth := 1
                linecolor := Colors.RGB(0.75*[1,1,1]...)
                Niϕϕ(k,T2),Ni_Gaϕ(k,T2)
            end
            Offset=k-180*floor(Int,k/180);
            if sign(sind(k))==1
                TextX=Niϕϕ(k,Ni_La(180-Offset))
                TextY=Ni_Gaϕ(k,Ni_La(180-Offset))
            else
                TextX=Niϕϕ(k,-Ni_La(Offset))+360
                TextY=Ni_Gaϕ(k,-Ni_La(Offset))
            end
        end
        TextX
        annotations := (TextX,TextY,("$(string(k))°"))

        title --> "Nichols chart"
        grid --> false
        legend --> false
        xguide --> "Phase [deg]"
        yguide --> "Magnitude [dB]"

    end

    extremas = extrema(Gains)
    # colors = [:blue, :cyan, :green, :yellow, :orange, :red, :magenta]
    for (sysi,s) = enumerate(systems)
        ℜresp, ℑresp        = nyquist(s, w)[1:2]
        ℜdata               = dropdims(ℜresp, dims=(1,2))
        ℑdata               = dropdims(ℑresp, dims=(1,2))
        mag                 = 20*log10.(sqrt.(ℜdata.^2 + ℑdata.^2))
        angles              = 180/π*angle(im*ℑdata.+ℜdata)
        extremas = extrema([extremas..., extrema(mag)...])
        @series begin
            linewidth --> 2
            hover --> [Printf.@sprintf("ω = %.3f", w) for w in w]
            angles, mag
        end
    end
    ylims --> extremas
    xlims --> extrema(PCyc*360)
    nothing

end

@userplot Sigmaplot
"""
    sigmaplot(sys, args...; hz=false balance=true, extrema)
    sigmaplot(LTISystem[sys1, sys2...], args...; hz=false, balance=true, extrema)

Plot the singular values of the frequency response of the `LTISystem`(s). A
frequency vector `w` can be optionally provided.

- If `hz=true`, the plot x-axis will be displayed in Hertz, the input frequency vector is still treated as rad/s.
- `balance`: Call [`balance_statespace`](@ref) on the system before plotting.
- `extrema`: Only plot the largest and smallest singular values.

`kwargs` is sent as argument to Plots.plot.
"""
sigmaplot
@recipe function sigmaplot(p::Sigmaplot; hz=false, balance=true, extrema=false)
    systems, w = _processfreqplot(Val{:sigma}(), p.args...)
    ws = (hz ? 1/(2π) : 1) .* w
    ny, nu = size(systems[1])
    nw = length(w)
    title --> "Sigma Plot"
    xguide --> (hz ? "Frequency [Hz]" : "Frequency [rad/s]")
    yguide --> "Singular Values $_PlotScaleStr"
    for (si, s) in enumerate(systems)
        if balance
            s = balance_statespace(s)[1]
        end
        sv = sigma(s, w)[1]'
        if extrema && size(sv, 2) > 2
            sv = sv[:, [1, end]]
        end
        if _PlotScale == "dB"
            sv = 20*log10.(sv)
        end
        @series begin
            xscale --> :log10
            yscale --> _PlotScaleFunc
            seriescolor --> si
            _to1series(ws, sv)
        end
    end
end

"This is a helper function to make multiple series into one series separated by `Inf`. This makes plotting vastly more efficient. It's also useful to make many lines appear as a single series with a single legend entry."
function _to1series(x,y)
    r,c = size(y)
    y2 = vec([y; fill(Inf, 1, c)])
    x2 = repeat([x; Inf], c)
    x2,y2
end

_to1series(y) = _to1series(1:size(y,3),y)

@userplot Marginplot
"""
    fig = marginplot(sys::LTISystem [,w::AbstractVector];  balance=true, kwargs...)
    marginplot(sys::Vector{LTISystem}, w::AbstractVector;  balance=true, kwargs...)

Plot all the amplitude and phase margins of the system(s) `sys`.

- A frequency vector `w` can be optionally provided.
- `balance`: Call [`balance_statespace`](@ref) on the system before plotting.
- `adjust_phase_start`: If true, the phase will be adjusted so that it starts at -90*intexcess degrees, where `intexcess` is the integrator excess of the system.

`kwargs` is sent as argument to RecipesBase.plot.
"""
marginplot
@recipe function marginplot(p::Marginplot; plotphase=true, hz=false, balance=true, adjust_phase_start=true, adaptive=true)
    systems, w = _processfreqplot(Val{:bode}(), p.args...; adaptive)
    ny, nu = size(systems[1])
    s2i(i,j) = LinearIndices((nu,(plotphase ? 2 : 1)*ny))[j,i]
    layout --> ((plotphase ? 2 : 1)*ny, nu)
    titles = Array{AbstractString}(undef, nu,ny,2,2)
    titles[:,:,1,1] .= _PlotScale == "dB" ? "Gm (dB): " : "Gm: "
    titles[:,:,2,1] .= "Pm: "
    titles[:,:,1,2] .= "ω gm: "
    titles[:,:,2,2] .= "ω pm: "
    layout --> (2ny, nu)
    label --> ""
    for (si, s) in enumerate(systems)
        if balance
            s = balance_statespace(s)[1]
        end
        bmag, bphase = bode(s, w)

        for j=1:nu
            for i=1:ny
                wgm, gm, wpm, pm, fullPhase, phasedata = sisomargin(s[i,j],w; full=true, allMargins=true, adjust_phase_start)
                if length(gm) > 5
                    @warn "Only showing smallest 5 out of $(length(gm)) gain margins"
                    idx = sortperm(gm)
                    wgm = wgm[idx[1:5]]
                    gm = gm[idx[1:5]]
                end
                # Let's be reasonable, only plot 5 smallest phase margins
                if length(pm) > 5
                    @warn "Only showing \"smallest\" 5 out of $(length(pm)) phase margins"
                    idx = sortperm(pm)
                    wgm = wpm[idx[1:5]]
                    gm = pm[idx[1:5]]
                end
                if _PlotScale == "dB"
                    mag = 20 .* log10.(1 ./ gm)
                    oneLine = 0
                    @. bmag = 20*log10(bmag)
                    titles[j,i,1,1] *= "["*join([Printf.@sprintf("%3.2g",20log10(v)) for v in gm],", ")*"] "
                else
                    mag = 1 ./ gm
                    oneLine = 1
                    titles[j,i,1,1] *= "["*join([Printf.@sprintf("%3.2g",v) for v in gm],", ")*"] "
                end
                titles[j,i,1,2] *= "["*join([Printf.@sprintf("%3.2g",v) for v in wgm],", ")*"] "


                subplot := min(s2i((plotphase ? (2i-1) : i),j), prod(plotattributes[:layout]))
                if si == length(systems)
                    title := (titles[j,i,1,1]*" "*titles[j,i,1,2])
                end
                @series begin
                    lab = _get_plotlabel(s, i, j)
                    if lab !== nothing
                        label --> lab
                    end
                    primary := true
                    seriestype := :bodemag
                    m = bmag[i, j, :]
                    if adaptive
                        lmag = _PlotScale == "dB" ? m : log.(m)
                        wsi, _, inds = downsample(w, lmag, _span(lmag)/500)
                        wsi, m[inds]
                    else
                        w, m
                    end
                end
                
                #Plot gain margins
                @series begin
                    primary := false
                    color --> :gray
                    linestyle --> :dash
                    [w[1],w[end]], [oneLine,oneLine]
                end
                @series begin
                    primary := false
                    [wgm wgm]', [ones(length(mag)) mag]'
                end
                plotphase || continue


                titles[j,i,2,1] *=  "["*join([Printf.@sprintf("%3.2g°",v) for v in pm],", ")*"] "
                titles[j,i,2,2] *=  "["*join([Printf.@sprintf("%3.2g",v) for v in wpm],", ")*"] "
                
                # Phase margins
                subplot := s2i(2i,j)
                if si == length(systems)
                    title := (titles[j,i,2,1]*" "*titles[j,i,2,2])
                end
                @series begin
                    primary := true
                    seriestype := :bodephase
                    if adaptive
                        downsample(w, phasedata, _span(phasedata)/500)[1:2]
                    else
                        w, phasedata
                    end
                end
                @series begin
                    primary := false
                    color --> :gray
                    linestyle --> :dash
                    seriestype := :hline
                    ((fullPhase .- pm) .* ones(1, 2))'
                end
                @series begin
                    primary := false
                    [wpm wpm]', [fullPhase fullPhase-pm]'
                end
            end
        end
    end
end

# HELPERS:

function _same_io_dims(systems::LTISystem...)
    sizes = map(size, systems)
    return all(s -> s == sizes[1], sizes)
end

function _default_time_data(systems::Vector{T}) where T<:LTISystem
    sample_times = [_default_dt(i) for i in systems]
    tfinal = 100*maximum(sample_times)
    return sample_times, tfinal
end
_default_time_data(sys::LTISystem) = _default_time_data(LTISystem[sys])


@userplot Pzmap
"""
    fig = pzmap(fig, system, args...; hz = false, kwargs...)

Create a pole-zero map of the `LTISystem`(s) in figure `fig`, `args` and `kwargs` will be sent to the `scatter` plot command.

To customize the unit-circle drawn for discrete systems, modify the line attributes, e.g., `linecolor=:red`.

If `hz` is true, all poles and zeros are scaled by 1/2π.
"""
pzmap
@recipe function pzmap(p::Pzmap; hz=false)
    systems = p.args[1]
    seriestype := :scatter
    framestyle --> :zerolines
    title --> "Pole-zero map"
    legend --> false
    for (i, system) in enumerate(systems)
        p = poles(system)
        z = tzeros(system)
        if hz
            p ./= 2π
            z ./= 2π
        end
        if !isempty(z)
            @series begin
                group --> i
                markershape --> :c
                markersize --> 15.
                markeralpha --> 0.5
                real(z),imag(z)
            end
        end
        if !isempty(p)
            @series begin
                group --> i
                markershape --> :xcross
                markersize --> 15.
                markeralpha --> 0.5
                real(p),imag(p)
            end
        end

        if isdiscrete(system)
            plots_id = Base.PkgId(UUID("91a5bcdd-55d7-5caf-9e0b-520d859cae80"), "Plots")
            haskey(Base.loaded_modules, plots_id) || error("Call using Plots before calling this function")
            Plots = Base.loaded_modules[plots_id]
            @series begin
                seriestype := :shape
                fillalpha := 0
                Plots.partialcircle(0, 2π, 100)
            end
        end
    end
end
pzmap(sys::LTISystem; kwargs...) = pzmap([sys]; kwargs...)
pzmap!(sys::LTISystem; kwargs...) = pzmap!([sys]; kwargs...)

"""
    fig = gangoffourplot(P::LTISystem, C::LTISystem; minimal=true, plotphase=false, Ms_lines = [1.0, 1.25, 1.5], Mt_lines = [], sigma = true, kwargs...)

Gang-of-Four plot.

`sigma` determines whether a [`sigmaplot`](@ref) is used instead of a [`bodeplot`](@ref) for MIMO `S` and `T`.
`kwargs` are sent as argument to RecipesBase.plot.
"""
function gangoffourplot(P::Union{<:Vector, LTISystem}, C::Vector, args...; minimal=true, Ms_lines = [1.0, 1.25, 1.5], Mt_lines = [], sigma = true,  plotphase=false, adaptive=true, kwargs...)    
    if P isa LTISystem # Don't broadcast over scalar (with size?)
        P = [P]
    end
    plots_id = Base.PkgId(UUID("91a5bcdd-55d7-5caf-9e0b-520d859cae80"), "Plots")
    haskey(Base.loaded_modules, plots_id) || error("Call using Plots before calling this function")
    Plots = Base.loaded_modules[plots_id]

    gofs = gangoffour.(P,C)
    S,D,N,T = ntuple(i->getindex.(gofs, i), 4)
    bp = (args...; kwargs...) -> sigma ? sigmaplot(args...; kwargs...) : bodeplot(args...; plotphase=false, adaptive, kwargs...)
    f1 = bp(S, args...; show=false, title="S = 1/(1+PC)", kwargs...)
    if !isnothing(Ms_lines) && !isempty(Ms_lines)
        Plots.hline!(Ms_lines', l=(:dash, [:green :orange :red :darkred :purple]), sp=1, primary=false, lab=string.(Ms_lines'), ylims=(1e-2,4))
    else
        Plots.hline!([1.0], l=(:dash, :black), sp=1, ylims=(1e-2,1.8))
    end
    f2 = bodeplot(D, args...; show=false, title="P/(1+PC)", plotphase=false, adaptive, kwargs...)
    Plots.hline!(ones(1, ninputs(D[1])*noutputs(D[1])), l=(:black, :dash), primary=false)
    f3 = bodeplot(N, args...; show=false, title="C/(1+PC)", plotphase=false, adaptive, kwargs...)
    f4 = bp(T, args...; show=false, title="T = PC/(1+PC)", ylims=(1e-2,4), kwargs...)
    if !isnothing(Mt_lines) && !isempty(Mt_lines)
        Plots.hline!(Mt_lines', l=(:dash, [:green :orange :red :darkred :purple]), primary=false, lab=string.(Mt_lines'), ylims=(1e-2,4))
    else
        Plots.hline!([1.0], l=(:dash, :black), ylims=(1e-2,4))
    end
    Plots.plot(f1,f2,f3,f4, ticks=:default, ylabel="", legend=:bottomright)
end


function gangoffourplot(P::Union{<:Vector, LTISystem},C::LTISystem, args...; kwargs...)
    gangoffourplot(P,[C], args...; kwargs...)
end

@userplot Rgaplot
"""
    rgaplot(sys, args...; hz=false)
    rgaplot(LTISystem[sys1, sys2...], args...; hz=false, balance=true)

Plot the relative-gain array entries of the `LTISystem`(s). A
frequency vector `w` can be optionally provided.

- If `hz=true`, the plot x-axis will be displayed in Hertz, the input frequency vector is still treated as rad/s.
- `balance`: Call [`balance_statespace`](@ref) on the system before plotting.

`kwargs` is sent as argument to Plots.plot.
"""
rgaplot
@recipe function rgaplot(p::Rgaplot; hz=false, balance=true)
    systems, w = _processfreqplot(Val{:sigma}(), p.args...)
    ws = (hz ? 1/(2π) : 1) .* w
    ny, nu = size(systems[1])
    nw = length(w)
    title --> "RGA Plot"
    xguide --> (hz ? "Frequency [Hz]" : "Frequency [rad/s]")
    yguide --> "Element magnitudes"
    for (si, s) in enumerate(systems)
        if balance
            s = balance_statespace(s)[1]
        end
        sv = abs.(relative_gain_array(s, w))
        for j in 1:size(sv, 1)
            for i in 1:size(sv, 2)
                @series begin
                    xscale --> :log10
                    label --> "System $si, from $i to $j"
                    ws, sv[j, i, :]
                end
            end
        end
    end
end


## Adaptive sampling
# Code adapted from https://github.com/iuliancioarca/AdaptiveSampling.jl/blob/master/LICENSE

function downsample(t,y,detail_th)
    # Compress signal by removing redundant points.
    # Adjust waveform detail/compression ratio with 'detail_th' (maximum allowed
    # difference between original and approximated points from the signal)
    yln           = length(y)
    idx_l         = 1
    idx_r         = yln
    idx_d_max     = 1
    cond_break    = true
    d_max         = 0.0
    M             = zeros(Int,yln) # hash table for relevant indices
    idx2save      = zeros(Int,yln+2)
    cnt           = 2
    idx2save[1:2] = [1,yln]
    while cond_break
        # get maximum error(difference) and index, between original chunk of signal
        # and linear approximation
        d_max, idx_d_max = get_d_max(idx_l,idx_r,y)
        # save all indices
        M[idx_d_max] = idx_r
        if d_max > detail_th
            # if computed error is greater than maximum allowed error, save
            # next point index and call get_d_max(idx_l,idx_r,y) at next
            # iteration; keep going towards leftmost branches
            cnt           = cnt + 1
            idx_r         = idx_d_max
            idx2save[cnt] = idx_d_max
        else
            # if computed error is smaller than maximum allowed error, stop, go
            # right(to the next waveform segment) and call get_d_max(idx_l,idx_r,y)
            # at the next iteration
            idx_l     = idx_r;
            if idx_l != yln
                idx_r = M[idx_l]
            else
                cond_break = false
            end
        end
    end
    # sort all indexes corresponding to relevent points and generate resampled
    # signal
    idx2save = idx2save[1:cnt]
    idx2save = sort(idx2save)
    t_new    = @view t[idx2save]
    y_new    = @view y[idx2save]
    return t_new, y_new, idx2save
end
function get_d_max(idx_l,idx_r,y)
    # cut segment to be resampled
    yp = view(y,idx_l:idx_r)
    # construct linear approximation
    dr = LinRange(y[idx_l], y[idx_r], length(yp))
    # compute distance(error) and get index of maximum error
    # -> this will be used for further splitting the
    # signal and will be part of the final resampled signal
    d_max     = 0.0
    idx_d_max = 1
    err_val   = 0.0
    for i = 1:length(yp)
        err_val = abs(yp[i] - dr[i])
        if err_val > d_max
            d_max     = err_val
            idx_d_max = i
        end
    end
    idx_d_max = idx_d_max + idx_l - 1
    return d_max, idx_d_max
end