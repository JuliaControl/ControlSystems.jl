import PyPlot
export lsimplot, stepplot, impulseplot, bodeplot, nyquistplot, sigmaplot, marginplot, setPlotScale

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
function lsimplot(systems::Vector{LTISystem}, u::AbstractVecOrMat,
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
        function $funcname(systems::Vector{LTISystem}, Ts_list::Vector, Tf::Real)
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
        $funcname(systems::Vector{LTISystem}, Tf::Real) =
                $funcname(systems, map(_default_Ts, systems), Tf)
        $funcname(systems::Vector{LTISystem}) =
                $funcname(systems, _default_time_data(systems)...)
        $funcname(systems::Vector{LTISystem}, t::AbstractVector) =
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
function bodeplot(systems::Vector{LTISystem}, w::AbstractVector)
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
                axes[2*i, j][:semilogx](w, phasedata)
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
bodeplot(systems::Vector{LTISystem}) =
    bodeplot(systems, _default_freq_vector(systems, :bode))
bodeplot(sys::LTISystem, args...) = bodeplot(LTISystem[sys], args...)

@doc """`nyquistplot(sys, args...)`, `nyquistplot(LTISystem[sys1, sys2...],
args...)`

Create a Nyquist plot of the `LTISystem`(s). A frequency vector `w` can be
optionally provided.""" ->
function nyquistplot(systems::Vector{LTISystem}, w::AbstractVector)
    if !_same_io_dims(systems...)
        error("All systems must have the same input/output dimensions")
    end
    ny, nu = size(systems[1])
    nw = length(w)
    fig, temp = PyPlot.subplots(ny, nu, sharex="col", sharey="row")
    # Ensure that `axes` is always a matrix of handles
    axes = ny == 1 ? reshape([temp], ny, nu) : temp
    for s = systems
        re_resp, im_resp = nyquist(s, w)[1:2]
        for j=1:nu
            for i=1:ny
                redata = reshape(re_resp[i, j, :], nw)
                imdata = reshape(im_resp[i, j, :], nw)
                line = axes[i, j][:plot](redata, imdata)[1]
                color = line[:get_color]()
                # Plot the mirror
                ax = axes[i, j]
                ax[:plot](redata, -imdata, color=color)
                # Add arrows at the midpoint
                mp = div(nw, 2)
                ax[:arrow](redata[mp], imdata[mp], redata[mp + 1] - redata[mp],
                        imdata[mp + 1] - imdata[mp], color=color, width=0.003)
                ax[:arrow](redata[mp], -imdata[mp], redata[mp - 1] - redata[mp],
                        -imdata[mp - 1] + imdata[mp], color=color, width=0.003)
            end
        end
    end
    # Add labels and titles
    fig[:suptitle]("Nyquist Plot", size=16)
    if ny*nu != 1
        for i=1:ny
            axes[i, 1][:set_ylabel]("To: y($i)", size=12, color="0.30")
        end
        for j=1:nu
            ax = axes[1, j]
            ax[:set_title]("From: u($j)", size=12, color="0.30")
            # Ensure the x axis includes -1
            xlims = ax[:get_xlim]()
            ax[:set_xlim]([min(-1, xlims[1]), xlims[2]])
        end
    end
    fig[:text](0.06, 0.5, "Imaginary Axis", ha="center", va="center",
            rotation="vertical", size=14)
    fig[:text](0.5, 0.04, "Real Axis", ha="center", va="center", size=14)
    # Add axis ticks
    for ax in axes
        ax[:set_yticks]([0.0], minor=true)
        ax[:yaxis][:grid](true, which="minor")
        ax[:set_xticks]([0.0], minor=true)
        ax[:xaxis][:grid](true, which="minor")
    end
    PyPlot.draw()
    return fig
end
nyquistplot(systems::Vector{LTISystem}) =
    nyquistplot(systems, _default_freq_vector(systems, :nyquist))
nyquistplot(sys::LTISystem, args...) = nyquistplot(LTISystem[sys], args...)

@doc """`sigmaplot(sys, args...)`, `sigmaplot(LTISystem[sys1, sys2...],
args...)`

Plot the singular values of the frequency response of the `LTISystem`(s). A
frequency vector `w` can be optionally provided.""" ->
function sigmaplot(systems::Vector{LTISystem}, w::AbstractVector)
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
sigmaplot(systems::Vector{LTISystem}) =
    sigmaplot(systems, _default_freq_vector(systems, :sigma))
sigmaplot(sys::LTISystem, args...) = sigmaplot(LTISystem[sys], args...)


function marginplot(systems::Vector{LTISystem}, w::AbstractVector)
    if !_same_io_dims(systems...)
        error("All systems must have the same input/output dimensions")
    end
    ny, nu = size(systems[1])
    fig = bodeplot(systems,w)
    ax = fig[:axes]
    
    for s = systems
        for j=1:nu
            for i=1:ny
                wgm, gm, wpm, pm, fullPhase = margin(s[i,j],w, full=true)
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
                for k=1:length(wpm)
                    ax[nu*(2*i-1)+j][:semilogx]([wpm[k];wpm[k]],[fullPhase[k];fullPhase[k]-pm[k]])
                    ax[nu*(2*i-1)+j][:axhline](fullPhase[k]-pm[k],linestyle="--",color="gray")
                end
            end
        end
    end
    PyPlot.draw()
    return fig
end
marginplot(systems::Vector{LTISystem}) =
    marginplot(systems, _default_freq_vector(systems, :bode))
marginplot(sys::LTISystem, args...) = marginplot(LTISystem[sys], args...)


# HELPERS:

function _same_io_dims(systems::LTISystem...)
    sizes = map(size, systems)
    return all(s -> s == sizes[1], sizes)
end

function _default_time_data(systems::Vector{LTISystem})
    sample_times = [_default_Ts(i) for i in systems]
    Tf = 100*maximum(sample_times)
    return sample_times, Tf
end
_default_time_data(sys::LTISystem) = _default_time_data(LTISystem[sys])
