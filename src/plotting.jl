import PyPlot
export lsimplot, stepplot, impulseplot

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
    fig, axes = PyPlot.subplots(ny, 1) 
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
            fig, axes = PyPlot.subplots(ny, nu) 
            if ny*nu == 1
                axes = [axes]
            end
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


# HELPERS:

function _same_io_dims(systems::LTISystem...)
    sizes = map(size, systems)
    return all([s == sizes[1] for s in sizes])
end

function _default_time_data(systems::Vector{LTISystem})
    sample_times = [_default_Ts(i) for i in systems]
    Tf = 100*maximum(sample_times)
    return sample_times, Tf
end
_default_time_data(sys::LTISystem) = _default_time_data(LTISystem[sys])
