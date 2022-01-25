import Colors
export bodeplot, nyquistplot, sigmaplot, marginplot, setPlotScale, gangoffour, gangoffourplot, gangofseven, pzmap, pzmap!, nicholsplot

_PlotScale = "log10"
_PlotScaleFunc = :log10
_PlotScaleStr = ""



"""
    setPlotScale(str)

Set the default scale of magnitude in `bodeplot` and `sigmaplot`.
`str` should be either `"dB"` or `"log10"`.
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
# Get atributes from xlims or ylims
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
#             minorText = [" " for t in minorText]#fill!(minorText, L" ")
#         end
#         perm = sortperm([major; minor])
#         return [major; minor][perm], [majorText; minorText][perm]

#     else
#         return major, majorText
#     end
# end


# This will be called on plot(lsim(sys, args...))
@recipe function simresultplot(r::SimResult; plotu=false, plotx=false)
    ny, nu, nx = r.ny, r.nu, r.nx
    t = r.t
    n_series = size(r.y, 3) # step and impulse produce multiple results
    nplots = ny
    plotu && (nplots += nu)
    plotx && (nplots += nx)
    layout --> (nplots, 1)
    seriestype --> (iscontinuous(r.sys) ? :path : :steppost)
    for ms in 1:n_series
        for i=1:ny
            ytext = (ny > 1) ? "y($i)" : "y"
            @series begin
                xguide  --> "Time (s)"
                yguide  --> ytext
                label   --> (n_series > 1 ? "From u($(ms))" : "")
                subplot --> i
                t,  r.y[i, :, ms]
            end
        end
    end 
    plotind = ny+1
    if plotu # bug in recipe system, can't use `plotu || return`
        for i=1:nu
            utext = (nu > 1) ? "u($i)" : "u"
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
            xtext = (nx > 1) ? "x($i)" : "x"
            @series begin
                xguide  --> "Time (s)"
                yguide  --> xtext
                subplot --> plotind
                label --> ""
                t,  r.x[i, :]
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

"""
    _processfreqplot(plottype, system::LTISystem, [w])
    _processfreqplot(plottype, system::AbstractVector{<:LTISystem}, [w])

Calculate default frequency vector and put system in array of not already array.
`plottype` is one of `Val{:bode}, Val{:nyquist}, ...`
for which `_default_freq_vector` is defined.
Check that system dimensions are compatible.
"""
_processfreqplot(plottype, system::LTISystem, args...) =
    _processfreqplot(plottype, [system], args...)
# Catch when system is not vector, with and without frequency input

# Cantch correct form
function _processfreqplot(plottype, systems::AbstractVector{<:LTISystem},
            w = _default_freq_vector(systems, plottype))

    if !_same_io_dims(systems...)
        error("All systems must have the same input/output dimensions")
    end
    return systems, w
end


@userplot Bodeplot
## FREQUENCY PLOTS ##
"""
    fig = bodeplot(sys, args...)
    bodeplot(LTISystem[sys1, sys2...], args...; plotphase=true, kwargs...)

Create a Bode plot of the `LTISystem`(s). A frequency vector `w` can be
optionally provided. To change the Magnitude scale see `setPlotScale(str)`
                                            
If `hz=true`, the plot x-axis will be displayed in Hertz, the input frequency vector is still treated as rad/s.

`kwargs` is sent as argument to RecipesBase.plot.
"""
bodeplot

@recipe function bodeplot(p::Bodeplot; plotphase=true, ylimsphase=(), unwrap=true, hz=false)
    systems, w = _processfreqplot(Val{:bode}(), p.args...)
    ws = (hz ? 1/(2π) : 1) .* w
    ny, nu = size(systems[1])
    s2i(i,j) = LinearIndices((nu,(plotphase ? 2 : 1)*ny))[j,i]
    layout --> ((plotphase ? 2 : 1)*ny, nu)
    nw = length(w)
    # xticks --> getLogTicks(ws, getlims(:xlims, plotattributes, ws))
    grid   --> true

    for (si,s) = enumerate(systems)
        mag, phase = bode(s, w)[1:2]
        if _PlotScale == "dB" # Set by setPlotScale(str) globally
            mag = 20*log10.(mag)
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
                @series begin
                    yscale    --> _PlotScaleFunc
                    xscale    --> :log10
                    # if _PlotScale != "dB"
                    #     yticks    --> getLogTicks(magdata, getlims(:ylims, plotattributes, magdata))
                    # end
                    xguide    --> xlab
                    yguide    --> "Magnitude $_PlotScaleStr"
                    subplot   --> min(s2i((plotphase ? (2i-1) : i),j), prod(plotattributes[:layout]))
                    # title     --> "Bode plot from: u($j)"
                    # label     --> "\$G_{$(si)}\$"
                    group     --> group_ind
                    ws, magdata
                end
                plotphase || continue

                @series begin
                    xscale    --> :log10
                    # ylims      := ylimsphase
                    # yticks    --> yphaseticks
                    yguide    --> "Phase (deg)"
                    subplot   --> s2i(2i,j)
                    xguide    --> (hz ? "Frequency [Hz]" : "Frequency [rad/s]")
                    # label     --> "\$G_{$(si)}\$"
                    group     --> group_ind
                    ws, unwrap ? ControlSystems.unwrap(phasedata.*(pi/180)).*(180/pi) : phasedata
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
    yscale --> :log10
    xscale --> :log10
    yguide --> "Magnitude"
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
    fig = nyquistplot(sys; Ms_circles=Float64[], Mt_circles=Float64[], unit_circle=false, hz=false, critical_point=-1, kwargs...)
    nyquistplot(LTISystem[sys1, sys2...]; Ms_circles=Float64[], Mt_circles=Float64[], unit_circle=false, hz=false, critical_point=-1, kwargs...)

Create a Nyquist plot of the `LTISystem`(s). A frequency vector `w` can be
optionally provided.

- `unit_circle`: if the unit circle should be displayed. The Nyquist curve crosses the unit circle at the gain corssover frequency.
- `Ms_circles`: draw circles corresponding to given levels of sensitivity (circles around -1 with  radii `1/Ms`). `Ms_circles` can be supplied as a number or a vector of numbers. A design staying outside such a circle has a phase margin of at least `2asin(1/(2Ms))` rad and a gain margin of at least `Ms/(Ms-1)`.
- `Mt_circles`: draw circles corresponding to given levels of complementary sensitivity. `Mt_circles` can be supplied as a number or a vector of numbers.
- `critical_point`: point on real axis to mark as critical for encirclements

If `hz=true`, the hover information will be displayed in Hertz, the input frequency vector is still treated as rad/s.

`kwargs` is sent as argument to plot.
"""
nyquistplot
@recipe function nyquistplot(p::Nyquistplot; Ms_circles=Float64[], Mt_circles=Float64[], unit_circle=false, hz=false, critical_point=-1)
    systems, w = _processfreqplot(Val{:nyquist}(), p.args...)
    ny, nu = size(systems[1])
    nw = length(w)
    layout --> (ny,nu)
    framestyle --> :zerolines
    s2i(i,j) = LinearIndices((nu,ny))[j,i]
    θ = range(0, stop=2π, length=100)
    S, C = sin.(θ), cos.(θ)
    for (si,s) = enumerate(systems)
        re_resp, im_resp = nyquist(s, w)[1:2]
        for j=1:nu
            for i=1:ny
                redata = re_resp[i, j, :]
                imdata = im_resp[i, j, :]
                ylims --> (min(max(-20,minimum(imdata)),-1), max(min(20,maximum(imdata)),1))
                xlims --> (min(max(-20,minimum(redata)),-1), max(min(20,maximum(redata)),1))
                title --> "Nyquist plot from: u($j)"
                yguide --> "To: y($i)"
                @series begin
                    subplot --> s2i(i,j)
                    label --> "\$G_{$(si)}\$"
                    hover --> [hz ? Printf.@sprintf("f = %.3f", w/2π) : Printf.@sprintf("ω = %.3f", w) for w in w]
                    (redata, imdata)
                end                
                
                if si == length(systems)
                    @series begin # Mark the critical point
                        subplot --> s2i(i,j)
                        primary := false
                        markershape := :xcross
                        seriescolor := :red
                        markersize := 5
                        seriesstyle := :scatter
                        [critical_point], [0]
                    end
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
                            (-1 .+ (1/Mt) * C, (1/Mt) * S)
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
    sigmaplot(sys, args...; hz=false)
    sigmaplot(LTISystem[sys1, sys2...], args...; hz=false)

Plot the singular values of the frequency response of the `LTISystem`(s). A
frequency vector `w` can be optionally provided.

If `hz=true`, the plot x-axis will be displayed in Hertz, the input frequency vector is still treated as rad/s.

`kwargs` is sent as argument to Plots.plot.
"""
sigmaplot
@recipe function sigmaplot(p::Sigmaplot; hz=false)
    systems, w = _processfreqplot(Val{:sigma}(), p.args...)
    ws = (hz ? 1/(2π) : 1) .* w
    ny, nu = size(systems[1])
    nw = length(w)
    title --> "Sigma Plot"
    xguide --> (hz ? "Frequency [Hz]" : "Frequency [rad/s]")
    yguide --> "Singular Values $_PlotScaleStr"
    for (si, s) in enumerate(systems)
        sv = sigma(s, w)[1]'
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

_to1series(y) = _to1series(1:size(y,1),y)

@userplot Marginplot
"""
    fig = marginplot(sys::LTISystem [,w::AbstractVector];  kwargs...)
    marginplot(sys::Vector{LTISystem}, w::AbstractVector;  kwargs...)

Plot all the amplitude and phase margins of the system(s) `sys`.
A frequency vector `w` can be optionally provided.

`kwargs` is sent as argument to RecipesBase.plot.
"""
@recipe function marginplot(p::Marginplot; plotphase=true)
    systems, w = _processfreqplot(Val{:bode}(), p.args...)
    ny, nu = size(systems[1])
    s2i(i,j) = LinearIndices((nu,(plotphase ? 2 : 1)*ny))[j,i]
    layout --> ((plotphase ? 2 : 1)*ny, nu)
    titles = Array{AbstractString}(undef, nu,ny,2,2)
    titles[:,:,1,1] .= "Gm: "
    titles[:,:,2,1] .= "Pm: "
    titles[:,:,1,2] .= "Wgm: "
    titles[:,:,2,2] .= "Wpm: "
    layout --> (2ny, nu)
    label --> ""
    for (si, s) in enumerate(systems)
        bmag, bphase = bode(s, w)
        for j=1:nu
            for i=1:ny
                wgm, gm, wpm, pm, fullPhase = sisomargin(s[i,j],w, full=true, allMargins=true)
                # Let's be reasonable, only plot 5 smallest gain margins
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
                else
                    mag = 1 ./ gm
                    oneLine = 1
                end
                titles[j,i,1,1] *= "["*join([Printf.@sprintf("%2.2f",v) for v in gm],", ")*"] "
                titles[j,i,1,2] *= "["*join([Printf.@sprintf("%2.2f",v) for v in wgm],", ")*"] "
                titles[j,i,2,1] *=  "["*join([Printf.@sprintf("%2.2f",v) for v in pm],", ")*"] "
                titles[j,i,2,2] *=  "["*join([Printf.@sprintf("%2.2f",v) for v in wpm],", ")*"] "
                
                @series begin
                    primary := true
                    subplot --> min(s2i((plotphase ? (2i-1) : i),j), prod(plotattributes[:layout]))
                    seriestype := :bodemag
                    w, bmag[i, j, :]
                end
                
                #Plot gain margins
                primary --> false
                @series begin
                    subplot --> min(s2i((plotphase ? (2i-1) : i),j), prod(plotattributes[:layout]))
                    color --> :gray
                    linestyle --> :dash
                    [w[1],w[end]], [oneLine,oneLine]
                end
                @series begin
                    subplot --> min(s2i((plotphase ? (2i-1) : i),j), prod(plotattributes[:layout]))
                    title --> (titles[j,i,1,1]*" "*titles[j,i,1,2])
                    [wgm wgm]', [ones(length(mag)) mag]'
                end
                plotphase || continue

                # Phase margins
                @series begin
                    primary := true
                    subplot --> s2i(2i,j)
                    seriestype := :bodephase
                    w, bphase[i, j, :]
                end
                @series begin
                    subplot --> s2i(2i,j)
                    color --> :gray
                    linestyle --> :dash
                    [w[1],w[end]], [oneLine,oneLine]
                end
                @series begin
                    subplot --> s2i(2i,j)
                    [wpm wpm]', [fullPhase fullPhase-pm]'
                end
                @series begin
                    title --> (titles[j,i,2,1]*" "*titles[j,i,2,2])
                    subplot --> s2i(2i,j)
                    color --> :gray
                    linestyle --> :dash
                    [w[1] w[end]]', ((fullPhase .- pm) .* ones(1, 2))'
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
    fig = pzmap(fig, system, args...; kwargs...)

Create a pole-zero map of the `LTISystem`(s) in figure `fig`, `args` and `kwargs` will be sent to the `scatter` plot command.
"""
pzmap
@recipe function pzmap(p::Pzmap)
    systems = p.args[1]
    seriestype := :scatter
    framestyle --> :zerolines
    title --> "Pole-zero map"
    legend --> false
    for (i, system) in enumerate(systems)
        p = poles(system)
        z = tzeros(system)
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
            θ = range(0,stop=2π,length=100)
            S,C = sin.(θ),cos.(θ)
            seriestype := :path
            @series begin
                linecolor := :black
                C,S
            end
        end
    end
end
pzmap(sys::LTISystem; kwargs...) = pzmap([sys]; kwargs...)
pzmap!(sys::LTISystem; kwargs...) = pzmap!([sys]; kwargs...)

"""
    fig = gangoffourplot(P::LTISystem, C::LTISystem; minimal=true, plotphase=false, Ms_lines = [1.0, 1.1, 1.2], sigma = true, kwargs...)

Gang-of-Four plot.

`sigma` determines whether a [`sigmaplot`](@ref) is used instead of a [`bodeplot`](@ref) for MIMO `S` and `T`.
`kwargs` are sent as argument to RecipesBase.plot.
"""
function gangoffourplot(P::Union{<:Vector, LTISystem}, C::Vector, args...; minimal=true, Ms_lines = [1.0, 1.1, 1.2], sigma = true,  plotphase=false, kwargs...)    
    if P isa LTISystem # Don't broadcast over scalar (with size?)
        P = [P]
    end
    plots_id = Base.PkgId(UUID("91a5bcdd-55d7-5caf-9e0b-520d859cae80"), "Plots")
    haskey(Base.loaded_modules, plots_id) || error("Call using Plots before calling this function")
    Plots = Base.loaded_modules[plots_id]

    gofs = gangoffour.(P,C)
    S,D,N,T = ntuple(i->getindex.(gofs, i), 4)
    bp = (args...; kwargs...) -> sigma ? sigmaplot(args...; kwargs...) : bodeplot(args...; plotphase=false, kwargs...)
    f1 = bp(S, args...; show=false, title="S = 1/(1+PC)", kwargs...)
    if !isnothing(Ms_lines) && !isempty(Ms_lines)
        Plots.hline!(Ms_lines', l=(:dash, [:green :orange :red :darkred :purple]), sp=1, primary=false, lab=string.(Ms_lines'), ylims=(3e-3,8))
    else
        Plots.hline!([1.0], l=(:dash, :black), sp=1, ylims=(3e-3,8))
    end
    f2 = bodeplot(D, args...; show=false, title="P/(1+PC)", plotphase=false, kwargs...)
    Plots.hline!(ones(1, ninputs(D[1])*noutputs(D[1])), l=(:black, :dash), primary=false)
    f3 = bodeplot(N, args...; show=false, title="C/(1+PC)", plotphase=false, kwargs...)
    f4 = bp(T, args...; show=false, title="T = PC/(1+PC)", ylims=(3e-3,8), kwargs...)
    Plots.hline!([1], l=(:black, :dash), primary=false)
    Plots.plot(f1,f2,f3,f4, ticks=:default, ylabel="", legend=:bottomright)
end


function gangoffourplot(P::Union{<:Vector, LTISystem},C::LTISystem, args...; kwargs...)
    gangoffourplot(P,[C], args...; kwargs...)
end

@userplot Rgaplot
"""
    rgaplot(sys, args...; hz=false)
    rgaplot(LTISystem[sys1, sys2...], args...; hz=false)

Plot the relative-gain array entries of the `LTISystem`(s). A
frequency vector `w` can be optionally provided.

If `hz=true`, the plot x-axis will be displayed in Hertz, the input frequency vector is still treated as rad/s.

`kwargs` is sent as argument to Plots.plot.
"""
rgaplot
@recipe function rgaplot(p::Rgaplot; hz=false)
    systems, w = _processfreqplot(Val{:sigma}(), p.args...)
    ws = (hz ? 1/(2π) : 1) .* w
    ny, nu = size(systems[1])
    nw = length(w)
    title --> "RGA Plot"
    xguide --> (hz ? "Frequency [Hz]" : "Frequency [rad/s]")
    yguide --> "Element magnitudes"
    for (si, s) in enumerate(systems)
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
