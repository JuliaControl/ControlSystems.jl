import Colors
export lsimplot, stepplot, impulseplot, bodeplot, nyquistplot, sigmaplot, marginplot, setPlotScale, gangoffour, gangoffourplot, gangofseven, pzmap, pzmap!, nicholsplot

getColorSys(i,Nsys)   = convert(Colors.RGB,Colors.HSV(360*((i-1)/Nsys)^1.5,0.9,0.8))
function getStyleSys(i,Nsys)
    Ncmax = 5
    if Nsys <= Ncmax
        return Dict(:c => convert(Colors.RGB,Colors.HSV(360*((i-1)/Nsys)^1.5,0.9,0.8)) , :l => :solid)
    end
    styles = [:solid, :dash, :dot, :dashdot]
    Nstyles = min(ceil(Int,Nsys/Ncmax), length(styles))
    Nc = ceil(Int,Nsys / Nstyles)
    istyle = min(floor(Int,i/Nc)+1 , length(styles))
    c = convert(Colors.RGB,Colors.HSV(360*((mod(i-1,Nc)/Nc))^1.5,0.9,0.8))

    return Dict(:c => c, :l => styles[istyle])
end

_PlotScale = "dB"
_PlotScaleFunc = :identity
_PlotScaleStr = "(dB)"

@doc """`setPlotScale(str)`

Set the default scale of magnitude in `bodeplot` and `sigmaplot`.
`str` should be either `"dB"` or `"log10"`.""" ->
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

@doc """`fig = lsimplot(sys::StateSpace, u, t[, x0, method]; kwargs...), lsimplot(sys::TransferFunction,u,t[,method]; kwargs...)`

`lsimplot(StateSpace[sys1, sys2...], u, t[, x0, method]; kwargs...), lsimplot(TransferFunction[sys1, sys2...], u, t[, method]; kwargs...)`

Calculate the time response of the `LTISystem`(s) to input `u`. If `x0` is
ommitted, a zero vector is used.

Continuous time systems are discretized before simulation. By default, the
method is chosen based on the smoothness of the input signal. Optionally, the
`method` parameter can be specified as either `:zoh` or `:foh`.

`kwargs` is sent as argument to Plots.plot.""" ->
function lsimplot{T<:StateSpace}(systems::Vector{T}, u::Union{AbstractVecOrMat,Function},
    t::AbstractVector, x0::VecOrMat=zeros(systems[1].nx, 1),
    method::Symbol=_issmooth(u) ? :foh : :zoh; kwargs...)
    if !_same_io_dims(systems...)
        error("All systems must have the same input/output dimensions")
    end
    ny, nu = size(systems[1])
    fig = Plots.subplot(n=ny, nr=ny)
    for (si, s) in enumerate(systems)
        y = lsim(s, u, t, x0, method)[1]
        for i=1:ny
            ydata = reshape(y[:, i], size(t, 1))
            style = iscontinuous(s) ? :path : :steppost
            ytext = (ny > 1) ? "Amplitude to: y($i)": "Amplitude"
            Plots.plot!(fig[i,1], t, ydata, l=style, xlabel="Time (s)", ylabel=ytext, title="System Response",  lab="\$G_\{$(si)\}\$"; getStyleSys(si,length(systems))..., kwargs...)
        end
    end
    return fig
end
lsimplot(sys::LTISystem, u::Union{AbstractVecOrMat,Function}, t::AbstractVector, args...; kwargs...) =
    lsimplot(StateSpace[sys], u, t, args...; kwargs...)
lsimplot{T<:LTISystem}(sys::Vector{T}, u::Union{AbstractVecOrMat,Function}, t::AbstractVector, args...; kwargs...) =
    lsimplot(StateSpace[s for s in sys], u, t, args...; kwargs...)

for (func, title) = ((:step, "Step Response"), (:impulse, "Impulse Response"))
    funcname = Symbol("$(func)plot")
    @eval begin
        function $funcname{T<:LTISystem}(systems::Vector{T}, Ts_list::Vector, Tf::Real; kwargs...)
            if !_same_io_dims(systems...)
                error("All systems must have the same input/output dimensions")
            end
            ny, nu = size(systems[1])
            fig = Plots.subplot(n = ny*nu, nr = ny)
            for (si,(s, Ts)) in enumerate(zip(systems, Ts_list))
                t = 0:Ts:Tf
                y = ($func)(s, t)[1]
                for i=1:ny
                    for j=1:nu
                        ydata = reshape(y[:, i, j], size(t, 1))
                        style = iscontinuous(s) ? :path : :steppost
                        ttext = (nu > 1 && i==1) ? $title*" from: u($j) " : $title
                        ytext = (ny > 1 && j==1) ? "Amplitude to: y($i)": "Amplitude"
                        Plots.plot!(fig[i,j], t, ydata, l=style, xlabel="Time (s)", ylabel=ytext, title=ttext, lab="\$G_\{$(si)\}\$"; getStyleSys(si,length(systems))..., kwargs...)
                    end
                end
            end
            return fig
        end
        $funcname{T<:LTISystem}(systems::Vector{T}, Tf::Real; kwargs...) =
            $funcname(systems, map(_default_Ts, systems), Tf; kwargs...)
        $funcname{T<:LTISystem}(systems::Vector{T}; kwargs...) =
            $funcname(systems, _default_time_data(systems)...; kwargs...)
        $funcname{T<:LTISystem}(systems::Vector{T}, t::AbstractVector; kwargs...) =
            $funcname(systems, repmat([t[2] - t[1]], length(systems)), t[end]; kwargs...)
        $funcname(sys::LTISystem, args...; kwargs...) = $funcname(LTISystem[sys], args...; kwargs...)
    end
end

@doc """`fig = stepplot(sys, args...)`, `stepplot(LTISystem[sys1, sys2...], args...)`

Plot the `step` response of the `LTISystem`(s). A final time `Tf` or a
time vector `t` can be optionally provided.""" -> stepplot

@doc """`fig = impulseplot(sys, args...)`, `impulseplot(LTISystem[sys1, sys2...], args...)`

Plot the `impulse` response of the `LTISystem`(s). A final time `Tf` or a
time vector `t` can be optionally provided.""" -> impulseplot


## FREQUENCY PLOTS ##
@doc """`fig = bodeplot(sys, args...)`, `bodeplot(LTISystem[sys1, sys2...], args...; plotphase=true, kwargs...)`

Create a Bode plot of the `LTISystem`(s). A frequency vector `w` can be
optionally provided.

`kwargs` is sent as argument to Plots.plot.""" ->
function bodeplot{T<:LTISystem}(systems::Vector{T}, w::AbstractVector; plotphase=true, kwargs...)
    if !_same_io_dims(systems...)
        error("All systems must have the same input/output dimensions")
    end
    ny, nu = size(systems[1])
    fig = Plots.subplot(n=(plotphase?2:1)*ny*nu, nc=nu)
    nw = length(w)
    for (si,s) = enumerate(systems)
        mag, phase = bode(s, w)[1:2]
        if _PlotScale == "dB"
            mag = 20*log10(mag)
        end
        xlab = plotphase ? "" : "Frequency (rad/s)"
        for j=1:nu
            for i=1:ny
                magdata = vec(mag[:, i, j])
                if all(magdata .== -Inf)
                    # 0 system, don't plot anything
                    continue
                end
                phasedata = vec(phase[:, i, j])
                Plots.plot!(fig[(plotphase?(2i-1):i),j], w, magdata, grid=true, yscale=_PlotScaleFunc, xscale=:log10, xlabel=xlab, title="Bode plot from: u($j)", ylabel="Magnitude $_PlotScaleStr", lab="\$G_\{$(si)\}\$"; getStyleSys(si,length(systems))..., kwargs...)
                plotphase && Plots.plot!(fig[2i,j], w, phasedata, grid=true, xscale=:log10, ylabel="Phase (deg)",xlabel="Frequency (rad/s)"; getStyleSys(si,length(systems))..., kwargs...)
            end
        end
    end

    return fig
end
bodeplot{T<:LTISystem}(systems::Vector{T}; plotphase=true, kwargs...) =
    bodeplot(systems, _default_freq_vector(systems, :bode); plotphase=plotphase, kwargs...)
bodeplot(sys::LTISystem, args...; plotphase=true, kwargs...) = bodeplot([sys], args...; plotphase=plotphase, kwargs...)

@doc """`fig = nyquistplot(sys; kwargs...)`, `nyquistplot(LTISystem[sys1, sys2...]; kwargs...)`

Create a Nyquist plot of the `LTISystem`(s). A frequency vector `w` can be
optionally provided.

`kwargs` is sent as argument to Plots.plot.""" ->
function nyquistplot{T<:LTISystem}(systems::Vector{T}, w::AbstractVector; neg=false, kwargs...)
    if !_same_io_dims(systems...)
        error("All systems must have the same input/output dimensions")
    end
    ny, nu = size(systems[1])
    nw = length(w)
    fig = Plots.subplot(n=ny*nu, nc= nu)
    # Ensure that `axes` is always a matrix of handles
    for (si,s) = enumerate(systems)
        re_resp, im_resp = nyquist(s, w)[1:2]
        for j=1:nu
            for i=1:ny
                redata = re_resp[:, i, j]
                imdata = im_resp[:, i, j]
                ylim = (min(max(-20,minimum(imdata)),-1), max(min(20,maximum(imdata)),1))
                xlim = (min(max(-20,minimum(redata)),-1), max(min(20,maximum(redata)),1))
                Plots.plot!(fig[i, j],redata, imdata, title="Nyquist plot from: u($j)", ylabel="To: y($i)", ylims=ylim, xlims=xlim; getStyleSys(si,length(systems))..., kwargs...)

                if si == length(systems)
                    v = linspace(0,2π,100)
                    S,C = sin(v),cos(v)
                    Plots.plot!(fig[i, j],C,S,l=:dash,c=:black, lab="")
                    Plots.plot!(fig[i, j],C-1,S,l=:dash,c=:red, grid=true, lab="")
                    # neg && Plots.plot!(fig[i, j],redata, -imdata, args...)
                end
            end
        end
    end
    return fig
end

nyquistplot{T<:LTISystem}(systems::Vector{T}; kwargs...) =
    nyquistplot(systems, _default_freq_vector(systems, :nyquist); kwargs...)
nyquistplot(sys::LTISystem, args...; kwargs...) = nyquistplot([sys], args...; kwargs...)


@doc """
fig = `nicholsplot{T<:LTISystem}(systems::Vector{T}, w::AbstractVector; kwargs...)`

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
used to specify colors, line styles etc. using regular Plots.jl syntax

This function is based on code subject to the two-clause BSD licence
Copyright 2011 Will Robertson
Copyright 2011 Philipp Allgeuer

""" ->
function nicholsplot{T<:LTISystem}(systems::Vector{T}, w::AbstractVector;
    text=true,
    Gains = [12, 6, 3, 1, 0.5, -0.5, -1, -3, -6, -10, -20, -40, -60],
    pInc = 30,
    sat = 0.4,
    val = 0.85,
    fontsize = 10,
    kwargs...)
    if !_same_io_dims(systems...)
        error("All systems must have the same input/output dimensions")
    end
    ny, nu = size(systems[1])

    if !iscontinuous(systems[1])
        w_nyquist = 2π/systems[1].Ts
        w = w[w.<= w_nyquist]
    end
    nw = length(w)

    # Gain circle functions
    angle(x)        = unwrap(atan2(imag(x),real(x)))
    RadM(m)         = abs(m/(m^2-1))
    CentreM(m)      = m^2/(1-m^2)
    Ny(mdb,t)       = CentreM(10^(mdb/20))+RadM(10^(mdb/20)).*(cosd(t)+im.*sind(t))
    Niϕ(mdb,t)      = rad2deg((angle(Ny(mdb,t))))
    Ni_Ga(mdb,t)    = 20.*log10(abs(Ny(mdb,t)))

    # Phase circle functions
    Radϕ(ϕ)         = 1./(2.*abs(sind(ϕ)))
    Nyℜ(ϕ,t)        = -0.5+Radϕ(ϕ).*cosd(t+mod(ϕ,180)-90)
    Nyℑ(ϕ,t)        = 1./(2.*tand(ϕ))+Radϕ(ϕ).*sind(t+mod(ϕ,180)-90)
    Niϕϕ(ϕ,t)       = rad2deg((angle(Nyℜ(ϕ,t)+im*Nyℑ(ϕ,t))))+360*floor(ϕ/360)
    Ni_Gaϕ(ϕ,t)     = 20.*log10(abs(Nyℜ(ϕ,t)+im*Nyℑ(ϕ,t)))
    Ni_La(ϕ)        = 0.090*10^(ϕ/60)
    getColor(mdb)   = convert(Colors.RGB,Colors.HSV(360*((mdb-minimum(Gains))/(maximum(Gains)-minimum(Gains)))^1.5,sat,val))

    fig             = Plots.plot()
    megaangles      = vcat(map(s -> 180/π*angle(squeeze(freqresp(s, w)[1],(2,3))), systems)...)
    filter!(x-> !isnan(x), megaangles)
    PCyc            = Set{Int}(floor(Int,megaangles/360))
    PCyc            = sort(collect(PCyc))

    #  Gain circles
    for k=Gains
        ϕVals   =Niϕ(k,0:0.1:360)
        GVals   =Ni_Ga(k,0:0.1:360)
        for l in PCyc
            Plots.plot!(fig,ϕVals+l*360,GVals,c=getColor(k), grid=false)
            if text
                offset  = (l+1)*360
                TextX   = Niϕ(k,210)+offset
                TextY   = Ni_Ga(k,210)
                Plots.plot!(fig,ann=(TextX,TextY,Plots.text("$(string(k)) dB",fontsize)))
            end
        end
    end

    #  Phase circles
    Phi=PCyc[1]*360:pInc:PCyc[end]*360
    T1=logspace(-4,log10(180),300)
    T2=[T1; 360-flipdim(T1,1)]

    for k=Phi
        if abs(sind(k))<1e-3
            Plots.plot!(fig,[k,k],[-110,25],c=Colors.RGB(0.75*[1, 1, 1]...))
            if cosd(5)>0
                TextX=k
                TextY=1
            else
                TextX=k
                TextY=-46.5
            end
        else
            Plots.plot!(fig,Niϕϕ(k,T2),Ni_Gaϕ(k,T2),c=Colors.RGB(0.75*[1,1,1]...))
            Offset=k-180*floor(Int,k/180);
            if sign(sind(k))==1
                TextX=Niϕϕ(k,Ni_La(180-Offset))
                TextY=Ni_Gaϕ(k,Ni_La(180-Offset))
            else
                TextX=Niϕϕ(k,-Ni_La(Offset))+360
                TextY=Ni_Gaϕ(k,-Ni_La(Offset))
            end
        end
        if text
            Plots.plot!(fig,ann=(TextX,TextY,Plots.text("$(string(k))°",fontsize)))
        end

        Plots.plot!(fig, title="Nichols chart", grid=false, legend=false)

    end
    dKwargs = Dict(kwargs)
    LW = "linewidth" ∈ keys(dKwargs) ? pop!(dKwargs,"linewidth") : 2

    # colors = [:blue, :cyan, :green, :yellow, :orange, :red, :magenta]
    for (sysi,s) = enumerate(systems)
        ℜresp, ℑresp        = nyquist(s, w)[1:2]
        ℜdata               = squeeze(ℜresp, (2,3))
        ℑdata               = squeeze(ℑresp, (2,3))
        mag                 = 20*log10(sqrt(ℜdata.^2 + ℑdata.^2))
        angles              = 180/π*angle(im*ℑdata.+ℜdata)
        Plots.plot!(fig,angles, mag; linewidth = LW, getStyleSys(sysi,length(systems))..., kwargs...)
    end

    return fig
end

nicholsplot{T<:LTISystem}(systems::Vector{T};kwargs...) =
    nicholsplot(systems, _default_freq_vector(systems, :nyquist);kwargs...)
nicholsplot(sys::LTISystem, args...; kwargs...) = nicholsplot([sys],args...; kwargs...)

@doc """`sigmaplot(sys, args...)`, `sigmaplot(LTISystem[sys1, sys2...], args...)`

Plot the singular values of the frequency response of the `LTISystem`(s). A
frequency vector `w` can be optionally provided.

`kwargs` is sent as argument to Plots.plot.""" ->
function sigmaplot{T<:LTISystem}(systems::Vector{T}, w::AbstractVector; kwargs...)
    if !_same_io_dims(systems...)
        error("All systems must have the same input/output dimensions")
    end
    ny, nu = size(systems[1])
    nw = length(w)
    fig = Plots.plot()
    for (si, s) in enumerate(systems)
        sv = sigma(s, w)[1]
        if _PlotScale == "dB"
            sv = 20*log10(sv)
        end
        for i in 1:size(sv, 2)
            Plots.plot!(fig, w, sv[:, i], xscale=:log10, yscale=_PlotScaleFunc; getStyleSys(si,length(systems))..., kwargs...)
        end
    end
    Plots.plot!(fig, title="Sigma Plot", xlabel="Frequency (rad/s)",
        ylabel="Singular Values $_PlotScaleStr")
    return fig
end
sigmaplot{T<:LTISystem}(systems::Vector{T}; kwargs...) =
    sigmaplot(systems, _default_freq_vector(systems, :sigma); kwargs...)
sigmaplot(sys::LTISystem, args...; kwargs...) = sigmaplot([sys], args...; kwargs...)

@doc """`fig = marginplot(sys::LTISystem [,w::AbstractVector];  kwargs...)`, `marginplot(sys::Vector{LTISystem}, w::AbstractVector;  kwargs...)`

Plot all the amplitude and phase margins of the system(s) `sys`.
A frequency vector `w` can be optionally provided.

`kwargs` is sent as argument to Plots.plot.""" ->
function marginplot{T<:LTISystem}(systems::Vector{T}, w::AbstractVector; kwargs...)
    if !_same_io_dims(systems...)
        error("All systems must have the same input/output dimensions")
    end
    ny, nu = size(systems[1])
    fig = bodeplot(systems, w, kwargs...)

    titles = Array(AbstractString,nu,ny,2,2)
    titles[:,:,1,1] = "Gm: "
    titles[:,:,2,1] = "Pm: "
    titles[:,:,1,2] = "Wgm: "
    titles[:,:,2,2] = "Wpm: "
    for (si, s) in enumerate(systems)
        for j=1:nu
            for i=1:ny
                wgm, gm, wpm, pm, fullPhase = sisomargin(s[i,j],w, full=true, allMargins=true)
                if _PlotScale == "dB"
                    mag = 20*log10(1./gm)
                    oneLine = 0
                else
                    mag = 1./gm
                    oneLine = 1
                end
                for k=1:length(wgm)
                    #Plot gain margins
                    Plots.plot!(fig[2i-1,j], [wgm[k];wgm[k]], [1;mag[k]], lab=""; getStyleSys(si,length(systems))...)
                end
                #Plot gain line at 1
                Plots.plot!(fig[2i-1,j], [w[1],w[end]], [oneLine,oneLine], l=:dash, c=:gray, lab="")
                titles[j,i,1,1] *= "["*join([@sprintf("%2.2f",v) for v in gm],", ")*"] "
                titles[j,i,1,2] *= "["*join([@sprintf("%2.2f",v) for v in wgm],", ")*"] "
                for k=1:length(wpm)
                    #Plot the phase margins
                    Plots.plot!(fig[2i,j], [wpm[k];wpm[k]],[fullPhase[k];fullPhase[k]-pm[k]], lab=""; getStyleSys(si,length(systems))...)
                    #Plot the line at 360*k
                    Plots.plot!(fig[2i,j], [w[1],w[end]],(fullPhase[k]-pm[k])*ones(2), l=:dash, c=:gray, lab="")
                end
                titles[j,i,2,1] *=  "["*join([@sprintf("%2.2f",v) for v in pm],", ")*"] "
                titles[j,i,2,2] *=  "["*join([@sprintf("%2.2f",v) for v in wpm],", ")*"] "
            end
        end
    end
    for j = 1:nu
        for i = 1:ny
            Plots.title!(fig[2i-1,j], titles[j,i,1,1]*" "*titles[j,i,1,2])
            Plots.title!(fig[2i,j], titles[j,i,2,1]*" "*titles[j,i,2,2])
        end
    end
    return fig
end
marginplot{T<:LTISystem}(systems::Vector{T}; kwargs...) =
    marginplot(systems, _default_freq_vector(systems, :bode); kwargs...)
marginplot(sys::LTISystem, args...; kwargs...) = marginplot([sys], args...; kwargs...)


# HELPERS:

function _same_io_dims(systems::LTISystem...)
    sizes = map(size, systems)
    return all(s -> s == sizes[1], sizes)
end

function _default_time_data{T<:LTISystem}(systems::Vector{T})
    sample_times = [_default_Ts(i) for i in systems]
    Tf = 100*maximum(sample_times)
    return sample_times, Tf
end
_default_time_data(sys::LTISystem) = _default_time_data(LTISystem[sys])


@doc """`fig = pzmap!(fig, system, args...; kwargs...)`

Create a pole-zero map of the `LTISystem`(s) in figure `fig`, `args` and `kwargs` will be sent to the `scatter` plot command.""" ->
function pzmap!(fig, system::LTISystem, args...; kwargs...)
    if system.nu + system.ny > 2
        warn("pzmap currently only supports SISO systems. Only transfer function from u₁ to y₁ will be shown")
    end

    z,p,k = zpkdata(system)
    !isempty(z[1]) && Plots.scatter!(fig,real(z[1]),imag(z[1]),m=:c,markersize=15., markeralpha=0.5, args...; kwargs...)
    !isempty(p[1]) && Plots.scatter!(fig,real(p[1]),imag(p[1]),m=:x,markersize=15., args...; kwargs...)
    Plots.title!("Pole-zero map")

    if system.Ts > 0
        v = linspace(0,2π,100)
        S,C = sin(v),cos(v)
        Plots.plot!(fig,C,S,l=:dash,c=:black, grid=true)
    end
    Plots.plot!(fig,legend=false)

    return fig
end

@doc """`fig = pzmap(system, args...; kwargs...)`

Create a pole-zero map of the `LTISystem`(s), `args` and `kwargs` will be sent to the `scatter` plot command.""" ->
pzmap(system::LTISystem, args...; kwargs...) = pzmap!(Plots.plot(), system, args...; kwargs...)

@doc """`fig = gangoffourplot(P::LTISystem, C::LTISystem)`, `gangoffourplot(P::Union{Vector, LTISystem}, C::Vector; plotphase=false)`

Gang-of-Four plot.

`kwargs` is sent as argument to Plots.plot.""" ->
function gangoffourplot(P::Union{Vector, LTISystem}, C::Vector, args...; plotphase=false, kwargs...)
    S,D,N,T = gangoffour(P,C)
    fig = bodeplot(LTISystem[[S[i] D[i]; N[i] T[i]] for i = 1:length(C)], args..., plotphase=plotphase; kwargs...)
    lower = plotphase ? 3 : 2
    Plots.plot!(fig[1,1],title="\$S = 1/(1+PC)\$")
    Plots.plot!(fig[1,2],title="\$D = P/(1+PC)\$")
    Plots.plot!(fig[lower,1],title="\$N = C/(1+PC)\$")
    Plots.plot!(fig[lower,2],title="\$T = PC/(1+PC\$)")
    return fig
end


function gangoffourplot(P::LTISystem,C::LTISystem, args...; plotphase=false, kwargs...)
    gangoffourplot(P,[C], args...; plotphase=plotphase, kwargs...)
end
