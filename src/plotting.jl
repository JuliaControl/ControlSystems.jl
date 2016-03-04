import PyPlot, Colors
export lsimplot, stepplot, impulseplot, bodeplot, nyquistplot, sigmaplot, marginplot, setPlotScale, gangoffour, gangoffourplot, gangofseven, pzmap, pzmap!, nicholsplot

_PlotScale = "dB"
_PlotScaleFunc = :semilogx
_PlotScaleStr = "(dB)"

@doc """`setPlotScale(str)`

Set the default scale of magnitude in `bodeplot` and `sigmaplot`.
`str` should be either `"dB"` or `"log10"`.""" ->
function setPlotScale(str::AbstractString)
    if str == "dB"
        plotSettings = (str, :semilogx, "(dB)")
    elseif str == "log10"
        plotSettings = (str, :loglog, "")
    else
        error("Scale must be set to either \"dB\" or \"log10\"")
    end
    global _PlotScale, _PlotScaleFunc, _PlotScaleStr
    _PlotScale, _PlotScaleFunc, _PlotScaleStr = plotSettings
end

@doc """`lsimplot(sys, u, t[, x0, method])`

`lsimplot(LTISystem[sys1, sys2...], u, t[, x0, method])`

Calculate the time response of the `LTISystem`(s) to input `u`. If `x0` is
ommitted, a zero vector is used.

Continuous time systems are discretized before simulation. By default, the
method is chosen based on the smoothness of the input signal. Optionally, the
`method` parameter can be specified as either `:zoh` or `:foh`.""" ->
function lsimplot{T<:LTISystem}(systems::Vector{T}, u::AbstractVecOrMat,
        t::AbstractVector, x0::VecOrMat=zeros(systems[1].nx, 1),
        method::Symbol=_issmooth(u) ? :foh : :zoh)
    if !_same_io_dims(systems...)
        error("All systems must have the same input/output dimensions")
    end
    ny, nu = size(systems[1])
    fig, axes = PyPlot.subplots(ny, 1, sharex=true)
    if ny == 1
        axes = [axes]
    end
    for s = systems
        y = lsim(s, u, t, x0, method)[1]
        for i=1:ny
            ax = axes[i]
            ydata = reshape(y[:, i], size(t, 1))
            if iscontinuous(s)
                ax[:plot](t, ydata)
            else
                ax[:step](t, ydata, where="post")
            end
        end
    end
    # Add labels and titles
    fig[:suptitle]("System Response", size=16)
    if ny != 1
        for i=1:ny
            axes[i, 1][:set_ylabel]("To: y($i)", size=12, color="0.30")
        end
    end
    fig[:text](0.5, 0.04, "Time (s)", ha="center", va="center", size=14)
    fig[:text](0.06, 0.5, "Amplitude", ha="center", va="center",
        rotation="vertical", size=14)
    PyPlot.draw()
    return fig
end
lsimplot(sys::LTISystem, u::AbstractVecOrMat, t::AbstractVector, args...) =
    lsimplot(LTISystem[sys], u, t, args...)


for (func, title) = ((:step, "Step Response"), (:impulse, "Impulse Response"))
    funcname = Symbol("$(func)plot")
    @eval begin
        function $funcname{T<:LTISystem}(systems::Vector{T}, Ts_list::Vector, Tf::Real)
            if !_same_io_dims(systems...)
                error("All systems must have the same input/output dimensions")
            end
            ny, nu = size(systems[1])
            fig, temp = PyPlot.subplots(ny, nu, sharex="col", sharey="row")
            # Ensure that `axes` is always a matrix of handles
            axes = ny == 1 ? reshape([temp], ny, nu) : temp
            for (s, Ts) in zip(systems, Ts_list)
                t = 0:Ts:Tf
                y = ($func)(s, t)[1]
                for i=1:ny
                    for j=1:nu
                        ax = axes[i, j]
                        ydata = reshape(y[:, i, j], size(t, 1))
                        if iscontinuous(s)
                            ax[:plot](t, ydata)
                        else
                            ax[:step](t, ydata, where="post")
                        end
                    end
                end
            end
            # Add labels and titles
            fig[:suptitle]($title, size=16)
            if ny*nu != 1
                for i=1:ny
                    axes[i, 1][:set_ylabel]("To: y($i)", size=12, color="0.30")
                end
                for j=1:nu
                    axes[1, j][:set_title]("From: u($j)", size=12, color="0.30")
                end
            end
            fig[:text](0.5, 0.04, "Time (s)", ha="center", va="center", size=14)
            fig[:text](0.06, 0.5, "Amplitude", ha="center", va="center",
                rotation="vertical", size=14)
            PyPlot.draw()
            return fig
        end
        $funcname{T<:LTISystem}(systems::Vector{T}, Tf::Real) =
            $funcname(systems, map(_default_Ts, systems), Tf)
        $funcname{T<:LTISystem}(systems::Vector{T}) =
            $funcname(systems, _default_time_data(systems)...)
        $funcname{T<:LTISystem}(systems::Vector{T}, t::AbstractVector) =
            $funcname(systems, repmat([t[2] - t[1]], length(systems)), t[end])
        $funcname(sys::LTISystem, args...) = $funcname(LTISystem[sys], args...)
    end
end

@doc """`stepplot(sys, args...)`, `stepplot(LTISystem[sys1, sys2...], args...)`

Plot the `step` response of the `LTISystem`(s). A final time `Tf` or a
time vector `t` can be optionally provided.""" -> stepplot

@doc """`impulseplot(sys, args...)`, `impulseplot(LTISystem[sys1, sys2...], args...)`

Plot the `impulse` response of the `LTISystem`(s). A final time `Tf` or a
time vector `t` can be optionally provided.""" -> impulseplot


## FREQUENCY PLOTS ##
@doc """`bodeplot(sys, args...)`, `bodeplot(LTISystem[sys1, sys2...], args...)`

Create a Bode plot of the `LTISystem`(s). A frequency vector `w` can be
optionally provided.""" ->
function bodeplot{T<:LTISystem}(systems::Vector{T}, w::AbstractVector)
    if !_same_io_dims(systems...)
        error("All systems must have the same input/output dimensions")
    end
    ny, nu = size(systems[1])
    fig, axes = PyPlot.subplots(2*ny, nu, sharex="col", sharey="row")
    nw = length(w)
    for s = systems
        mag, phase = bode(s, w)[1:2]
        if _PlotScale == "dB"
            mag = 20*log10(mag)
        end
        for j=1:nu
            for i=1:ny
                magdata = vec(mag[i, j, :])
                if all(magdata .== -Inf)
                    # 0 system, don't plot anything
                    continue
                end
                phasedata = vec(phase[i, j, :])
                axes[2*i - 1, j][_PlotScaleFunc](w, magdata)
                axes[2*i - 1, j][:grid](true, which="both")
                axes[2*i, j][:semilogx](w, phasedata)
                axes[2*i, j][:grid](true, which="both")
            end
        end
    end
    # Add labels and titles
    fig[:suptitle]("Bode Plot", size=16)
    if ny*nu != 1
        for i=1:2*ny
            div(i+1, 2)
            axes[i, 1][:set_ylabel]("To: y($(div(i + 1, 2)))",
                size=12, color="0.30")
        end
        for j=1:nu
            axes[1, j][:set_title]("From: u($j)", size=12, color="0.30")
        end
        fig[:text](0.06, 0.5, "Phase (deg), Magnitude $_PlotScaleStr", ha="center",
            va="center", rotation="vertical", size=14)
    else
        axes[1, 1][:set_ylabel]("Magnitude $_PlotScaleStr", size=14)
        axes[2, 1][:set_ylabel]("Phase (deg)", size=14)
    end
    fig[:text](0.5, 0.04, "Frequency (rad/s)", ha="center",
    va="center", size=14)
    PyPlot.draw()
    return fig
end
bodeplot{T<:LTISystem}(systems::Vector{T}) =
    bodeplot(systems, _default_freq_vector(systems, :bode))
bodeplot(sys::LTISystem, args...) = bodeplot(LTISystem[sys], args...)

@doc """ `nyquistplot(sys, args...)`, `nyquistplot(LTISystem[sys1, sys2...], args...)`

Create a Nyquist plot of the `LTISystem`(s). A frequency vector `w` can be
optionally provided.""" ->
function nyquistplot{T<:LTISystem}(systems::Vector{T}, w::AbstractVector, args...; neg=false)
    if !_same_io_dims(systems...)
        error("All systems must have the same input/output dimensions")
    end
    ny, nu = size(systems[1])
    nw = length(w)
    fig = Plots.subplot(n=ny*nu, nc= nu)
    # Ensure that `axes` is always a matrix of handles
    for s = systems
        re_resp, im_resp = nyquist(s, w)[1:2]
        for j=1:nu
            for i=1:ny
                redata = reshape(re_resp[i, j, :], nw)
                imdata = reshape(im_resp[i, j, :], nw)
                Plots.plot!(fig[i, j],redata, imdata, title="From: u($j)", ylabel="To: y($i)", args...)

                v = linspace(0,2π,100)
                S,C = sin(v),cos(v)
                Plots.plot!(fig[i, j],C,S,l=:dash,c=:black, grid=true)
                # neg && Plots.plot!(fig[i, j],redata, -imdata, args...)
            end
        end
    end
    return fig
end

nyquistplot{T<:LTISystem}(systems::Vector{T}; kwargs...) =
    nyquistplot(systems, _default_freq_vector(systems, :nyquist); kwargs...)
nyquistplot(sys::LTISystem, args...; kwargs...) = nyquistplot(LTISystem[sys], args...; kwargs...)


@doc """
`nicholsplot{T<:LTISystem}(systems::Vector{T}, w::AbstractVector)`

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
    megaangles      = collect(map(s -> 180/π*angle(squeeze(freqresp(s, w)[1],(1,2))), systems)...)
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
    getColorSys(i)   = convert(Colors.RGB,Colors.HSV(360*((i-1)/(length(systems)))^1.5,0.9,0.8))
    for (sysi,s) = enumerate(systems)
        ℜresp, ℑresp        = nyquist(s, w)[1:2]
        ℜdata               = squeeze(ℜresp, (1,2))
        ℑdata               = squeeze(ℑresp, (1,2))
        mag                 = 20*log10(sqrt(ℜdata.^2 + ℑdata.^2))
        angles              = 180/π*angle(im*ℑdata.+ℜdata)
        Plots.plot!(fig,angles, mag; linewidth = LW, c = getColorSys(sysi), kwargs...)
    end

    return fig
end

# function feedback{T<:LTISystem}(systems::Vector{T})
#     systems2 = deepcopy(systems)
#     for (i,G) in enumerate(systems)
#         if isa(G,TransferFunction) #|| isa(G,zpkdata)
#             systems2[i] = minreal(G/(1+G))
#         elseif isa(G,StateSpace)
#             systems2[i].A = G.A-G.B*
#         end
#     end
# end

nicholsplot{T<:LTISystem}(systems::Vector{T};kwargs...) =
    nicholsplot(systems, _default_freq_vector(systems, :nyquist);kwargs...)
nicholsplot(sys::LTISystem, args...; kwargs...) = nicholsplot(LTISystem[sys],args...; kwargs...)


function unwrap(ϕ)
    if length(ϕ) == 1
        return ϕ
    end
    δϕ = diff(ϕ)
    δϕs = mod(δϕ+π,2π) - π
    δϕs[(δϕs.==-π) & (δϕ.>0)] = π
    δϕadd = δϕs - δϕ
    δϕadd[abs(δϕ).<π] = 0

    ϕ[2:end] = ϕ[2:end] + cumsum(δϕadd)
    return ϕ
end

@doc """`sigmaplot(sys, args...)`, `sigmaplot(LTISystem[sys1, sys2...], args...)`

Plot the singular values of the frequency response of the `LTISystem`(s). A
frequency vector `w` can be optionally provided.""" ->
function sigmaplot{T<:LTISystem}(systems::Vector{T}, w::AbstractVector)
    if !_same_io_dims(systems...)
        error("All systems must have the same input/output dimensions")
    end
    ny, nu = size(systems[1])
    nw = length(w)
    fig, ax = PyPlot.subplots(1, 1)
    for s = systems
        sv = sigma(s, w)[1]
        if _PlotScale == "dB"
            sv = 20*log10(sv)
        end
        # Plot the first singular value, grab the line color, then plot the
        # remaining values all in the same color.
        line = ax[_PlotScaleFunc](w, sv[1, :]')[1]
        color = line[:get_color]()
        for i in 2:size(sv, 1)
            ax[_PlotScaleFunc](w, sv[i, :]', color=color)
        end
    end
    ax[:set_title]("Sigma Plot", size=16)
    ax[:set_xlabel]("Frequency (rad/s)", size=14)
    ax[:set_ylabel]("Singular Values $_PlotScaleStr", size=14)
    PyPlot.draw()
    return fig
end
sigmaplot{T<:LTISystem}(systems::Vector{T}) =
    sigmaplot(systems, _default_freq_vector(systems, :sigma))
sigmaplot(sys::LTISystem, args...) = sigmaplot(LTISystem[sys], args...)


function marginplot{T<:LTISystem}(systems::Vector{T}, w::AbstractVector)
    if !_same_io_dims(systems...)
        error("All systems must have the same input/output dimensions")
    end
    ny, nu = size(systems[1])
    fig = bodeplot(systems,w)
    fig[:suptitle]("Margin Plot", size=16)

    ax = fig[:axes]
    titles = Array(String,nu,ny,2,2)
    titles[:,:,1,1] = "Gm: "
    titles[:,:,2,1] = "Pm: "
    titles[:,:,1,2] = "Wgm: "
    titles[:,:,2,2] = "Wpm: "
    for s = systems
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
                    ax[2*nu*(i-1)+j][_PlotScaleFunc]([wgm[k];wgm[k]],[1;mag[k]])
                end
                ax[2*nu*(i-1)+j][:axhline](oneLine,linestyle="--",color="gray")
                titles[nu,ny,1,1] *= "["*join([@sprintf("%2.2f",v) for v in gm],", ")*"] "
                titles[nu,ny,1,2] *= "["*join([@sprintf("%2.2f",v) for v in wgm],", ")*"] "
                for k=1:length(wpm)
                    ax[nu*(2*i-1)+j][:semilogx]([wpm[k];wpm[k]],[fullPhase[k];fullPhase[k]-pm[k]])
                    ax[nu*(2*i-1)+j][:axhline](fullPhase[k]-pm[k],linestyle="--",color="gray")
                end
                titles[nu,ny,2,1] *=  "["*join([@sprintf("%2.2f",v) for v in pm],", ")*"] "
                titles[nu,ny,2,2] *=  "["*join([@sprintf("%2.2f",v) for v in wpm],", ")*"] "
            end
        end
    end
    for j = 1:nu
        for i = 1:ny
            ax[2*nu*(i-1)+j][:set_title](titles[nu,ny,1,1]*" "*titles[nu,ny,1,2],loc="center")
            ax[nu*(2*i-1)+j][:set_title](titles[nu,ny,2,1]*" "*titles[nu,ny,2,2],loc="center")
            #if length(systems) > 1
            #        PyPlot.legend(["System "*string(k) for k in 1:length(systems)])
            #end
        end
    end
    PyPlot.draw()
    return fig
end
marginplot{T<:LTISystem}(systems::Vector{T}) =
    marginplot(systems, _default_freq_vector(systems, :bode))
marginplot(sys::LTISystem, args...) = marginplot(LTISystem[sys], args...)


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


# @doc """`pzmap(sys, args...)`, `pzmap(LTISystem[sys1, sys2...], args...)`
#
# Create a pole-zero map of the `LTISystem`(s).""" ->
# function pzmap(systems::Vector)
#     ny, nu = size(systems[1])
#     fig, axes = PyPlot.subplots(ny, nu, sharex="col", sharey="row")
#     for s = systems
#         z,p,k = zpkdata(s)
#         for j=1:nu
#             for i=1:ny
#                 axes[i - 1, j][plot](z,"o")
#                 axes[i - 1, j][plot](p,"x")
#             end
#         end
#     end
#     # Add labels and titles
#     fig[:suptitle]("Pole-zero map", size=16)
#     if ny*nu != 1
#         for i=1:ny
#             div(i+1, 2)
#             axes[i, 1][:set_ylabel]("To: y($(div(i + 1, 2)))",
#                     size=12, color="0.30")
#         end
#         for j=1:nu
#             axes[1, j][:set_title]("From: u($j)", size=12, color="0.30")
#         end
#     end
#     PyPlot.draw()
#     return fig
# end

@doc """`pzmap(sys)``

Create a pole-zero map of the `LTISystem`(s).""" ->
function pzmap!(fig, system::LTISystem, args...)
    if system.nu + system.ny > 2
        warn("pzmap currently only supports SISO systems. Only transfer function from u₁ to y₁ will be shown")
    end

    z,p,k = zpkdata(system)
    !isempty(z[1]) && Plots.scatter!(fig,real(z[1]),imag(z[1]),m=:c,markersize=15., markeralpha=0.5, args...)
    !isempty(p[1]) && Plots.scatter!(fig,real(p[1]),imag(p[1]),m=:x,markersize=15., args...)
    Plots.title!("Pole-zero map")

    if system.Ts > 0
        v = linspace(0,2π,100)
        S,C = sin(v),cos(v)
        Plots.plot!(fig,C,S,l=:dash,c=:black, grid=true)
    end
    Plots.plot!(fig,legend=false)

    return fig
end

pzmap(system::LTISystem, args...) = pzmap!(Plots.plot(), system::LTISystem, args...)

@doc """`gangoffourplot(P::LTISystem, C::LTISystem)`, `gangoffourplot(P::Union{Vector, LTISystem}, C::Vector)`

Gang-of-Four plot.""" ->
function gangoffourplot(P::Union{Vector, LTISystem}, C::Vector)
    S,D,N,T = gangoffour(P,C)
    fig = bodeplot(LTISystem[[S[i] D[i]; N[i] T[i]] for i = 1:length(C)])
    legend("S = \$1/(1+PC)\$","D = \$P/(1+PC)\$","N = \$C/(1+PC)\$","T = \$PC/(1+PC)\$")
    return fig
end

function gangoffourplot(P::LTISystem,C::LTISystem)
    S,D,N,T = gangoffour(P,C)
    fig = bodeplot([S D;N T])
    legend("S = \$1/(1+PC)\$","D = \$P/(1+PC)\$","N = \$C/(1+PC)\$","T = \$PC/(1+PC)\$")
    return fig
end
