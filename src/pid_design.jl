export pid, pidplots, rlocus, leadlink, laglink, leadlinkat, leadlinkcurve, stabregionPID, loopshapingPI

"""
Calculates and returns a PID controller on transfer function form.
`time` indicates whether or not the parameters are given as gains (default) or as time constants
`series` indicates  whether or not the series form or parallel form (default) is desired

`C = pid(; kp=0, ki=0; kd=0, time=false, series=false)`
"""
function pid(; kp=0, ki=0, kd=0, time=false, series=false)
    s = tf("s")
    if series
        return time ? kp*(1 + 1/(ki*s) + kd*s) : kp*(1 + ki/s + kd*s)
    else
        return time ? kp + 1/(ki*s) + kd*s : kp + ki/s + kd*s
    end
end

"""
Plots interesting figures related to closing the loop around process `P` with a PID controller
Send in a bunch of PID-parameters in any of the vectors kp, ki, kd. The vectors must be the same length.

`time` indicates whether or not the parameters are given as gains (default) or as time constants

`series` indicates  whether or not the series form or parallel form (default) is desired

Available plots are `:gof` for Gang of four, `:nyquist`, `:controller` for a bode plot of the controller TF and `:pz` for pole-zero maps

One can also supply a frequency vector ω to be used in Bode and Nyquist plots

`pidplots(P, args...; kps=0, kis=0, kds=0, time=false, series=false, ω=0)`

See also `loopshapingPI`, `stabregionPID`
"""
function pidplots(P::LTISystem, args...; kps=0, kis=0, kds=0, time=false, series=false, ω=0, grid = false)



    if grid
        kp = [i for i in kps, j in kis, k in kds][:]
        ki = [j for i in kps, j in kis, k in kds][:]
        kd = [k for i in kps, j in kis, k in kds][:]
        kps, kis, kds = kp, ki, kd
    else
        n = max(length(kps), length(kis), length(kds))
        kps = kps == 0 ? zeros(n) : kps
        kis = kis == 0 ? zeros(n) : kis
        kds = kds == 0 ? zeros(n) : kds
    end
    ω   = ω   == 0 ? logspace(-3,3,500) : ω

    getColorSys(i)   = convert(Colors.RGB,Colors.HSV(360*((i-1)/(length(kps)))^1.5,0.9,0.8))

    gof_        = in(:gof        ,args)
    nyquist_    = in(:nyquist    ,args)
    controller_ = in(:controller ,args)
    pz_         = in(:pz         ,args)
    nichols_    = in(:nichols    ,args)


    if nyquist_
        nq = Plots.plot()
    end
    if gof_
        bd = Plots.subplot(n=4,nc=2)
    end
    if pz_
        pz = Plots.plot()
    end
    if controller_
        cplot = plot()
    end

    for (i,kp) = enumerate(kps)
        ki = kis[i]
        kd = kds[i]
        label = "\$k_p = $(round(kp,3)), \\quad k_i = $(round(ki,3)), \\quad k_d = $(round(kd,3))\$"

        C = pid(kp=kp,ki=ki,kd=kd,time=time,series=series)
        S,D,N,T = gangoffour(P,C)

        if nyquist_
            NQ = nyquist(P*C,ω)
            redata = NQ[1][:]
            imdata = NQ[2][:]
            ylim = (max(-20,minimum(imdata)), min(20,maximum(imdata)))
            xlim = (max(-20,minimum(redata)), min(20,maximum(redata)))
            Plots.plot!(nq,redata,imdata, ylims=ylim, xlims=xlim, lab=label, c=getColorSys(i))
        end
        if gof_
            BD = bode(S,ω)
            Plots.plot!(bd[1,1],BD[3][:],BD[1][:], lab=label, c=getColorSys(i))
            BD = bode(D,ω)
            Plots.plot!(bd[1,2],BD[3][:],BD[1][:], lab=label, c=getColorSys(i))
            BD = bode(N,ω)
            Plots.plot!(bd[2,1],BD[3][:],BD[1][:], lab=label, c=getColorSys(i))
            BD = bode(T,ω)
            Plots.plot!(bd[2,2],BD[3][:],BD[1][:], lab=label, c=getColorSys(i))
        end
        if pz_
            pzmap!(pz,T)
        end
        if controller_
            BD = bode(C,ω)
            Plots.plot!(cplot,BD[3][:],BD[1][:], lab=label, c=getColorSys(i))
        end
    end

    nyquist_ && Plots.plot!(nq,legend=true, title="Nyquist curves")
    if gof_
        Plots.plot!(bd[1,1],legend=true, title="S", xscale=:log10, yscale=:log10)
        Plots.plot!(bd[1,2],legend=true, title="D", xscale=:log10, yscale=:log10)
        Plots.plot!(bd[2,1],legend=true, title="N", xscale=:log10, yscale=:log10)
        Plots.plot!(bd[2,2],legend=true, title="T", xscale=:log10, yscale=:log10)
    end
    if pz_
        Plots.plot!(pz,title="Pole-zero map")
    end
    if controller_
        Plots.plot!(cplot,title="Controller bode plot",legend=true, xscale=:log10, yscale=:log10)
    end



end

"""
`rlocus(P::LTISystem, K)` computes and plots the root locus of the SISO LTISystem P with a negative feedback loop and feedback gains `K`, if `K` is not provided, linspace(1e-6,500,10000) is used
"""
function rlocus(P::LTISystem, K=[])
    K = isempty(K) ? linspace(1e-6,500,10000) : K
    Z = tzero(P)
    poles = map(k -> pole(k*P/(1+k*P)), K)
    poles = cat(2,poles...)'
    redata = real(poles)
    imdata = imag(poles)
    ylim = (max(-50,minimum(imdata)), min(50,maximum(imdata)))
    xlim = (max(-50,minimum(redata)), min(50,maximum(redata)))
    Plots.plot(redata, imdata, legend=false,ylims=ylim, xlims=xlim)
    Plots.scatter!(real(Z), imag(Z), m=:c)
end

"""
`laglink(a, M; h=0)`

Returns a phase retarding link, the rule of thumb `a = 0.1ωc` guarantees less than 6 degrees phase margin loss. The bode curve will go from `M`, bend down at `a/M` and level out at 1 for frequencies > `a`
"""
function laglink(a, M; h=0)
    numerator = [1/a, 1]
    denominator = [M/a, 1]
    gain = M
    G = tf(gain*numerator,denominator)
    return  h <= 0 ? G : c2d(G,h)
end


"""
`leadlink(b, N, K; h=0)`

Returns a phase advancing link, the top of the phase curve is located at `ω = b√(N)` where the link amplification is `K√(N)` The bode curve will go from `K`, bend up at `b` and level out at `KN` for frequencies > `bN`

The phase advance at `ω = b√(N)` can be plotted as a function of `N` with `leadlinkcurve()`

Values of `N < 1` will give a phase retarding link.

See also `leadlinkat` `laglink`
"""
function leadlink(b, N, K; h=0)
    numerator = [1/b, 1]
    denominator = [1/(b*N), 1]
    gain = K
    G = tf(gain*numerator,denominator)
    return  h <= 0 ? G : c2d(G,h)

end

"""
`leadlinkat(ω, N, K; h=0)`

Returns a phase advancing link, the top of the phase curve is located at `ω` where the link amplification is `K√(N)` The bode curve will go from `K`, bend up at `ω/√(N)` and level out at `KN` for frequencies > `ω√(N)`

The phase advance at `ω` can be plotted as a function of `N` with `leadlinkcurve()`

Values of `N < 1` will give a phase retarding link.

See also `leadlink` `laglink`
"""
function leadlinkat(ω, N, K; h=0)
    b = ω / sqrt(N)
    return leadlink(b,N,K,h=h)
end

"""
Plot the phase advance as a function of `N` for a lead link (phase advance link)

If an input argument `s` is given, the curve is plotted from `s` to 10, else from 1 to 10.

See also `Leadlink, leadlinkat`
"""
function leadlinkcurve(start=1)
    N = linspace(start,10)
    dph = 180/pi*map(Ni->atan(sqrt(Ni))-atan(1/sqrt(Ni)), N)
    Plots.plot(N,dph, xlabel="N", ylabel="Phase advance [deg]")
end



"""
`fig, kp, ki = stabregionPID(P, [ω]; kd=0, doplot = true)`

Segments of the curve generated by this program
is the boundary of the stability region for a
process with transfer function P(s)
The PID controller is assumed to be on the form kp +ki/s +kd s

The curve is found by analyzing
P(s)\*C(s) = -1 ⟹\n
|PC| = |P| |C| = 1\n
arg(P) + arg(C) = -π


If `P` is a string (e.g. "exp(-sqrt(s))", the stability of feedback loops using PI-controllers can be analyzed for processes with models with arbitrary analytic functions

See also `stabregionPID`, `loopshapingPI`, `pidplots`
"""
function stabregionPID(P, ω = _default_freq_vector(P,:bode); kd=0, doplot = true)
    Pv      = squeeze(freqresp(P,ω)[1],(2,3))
    r       = abs(Pv)
    phi     = angle(Pv)
    kp      = -cos(phi)./r
    ki      = kd*ω.^2 - ω.*sin(phi)./r
    Plots.plot(kp,ki,linewidth = 1.5, xlabel="\$k_p\$", ylabel="\$k_i\$", title="Stability region of \$P, \\quad k_d = $(round(kd,4))\$"), kp, ki
end


function stabregionPID(P::AbstractString, ω = logspace(-3,1); kd=0, doplot = true)
    Pe      = parse(P)
    Pf(s)   = eval(:(s -> $(Pe)))(s)
    Pv      = Pf(im*ω)
    r       = abs(Pv)
    phi     = angle(Pv)
    kp      = -cos(phi)./r
    ki      = kd*ω.^2 - ω.*sin(phi)./r
    Plots.plot(kp,ki,linewidth = 1.5, xlabel="\$k_p\$", ylabel="\$k_i\$", title="Stability region of \$ $(replace(P,".","")), \\quad k_d = $(round(kd,4))\$"), kp, ki
end



"""
`kp,ki,C = loopshapingPI(P,ωp; ϕl,rl, phasemargin)`

Selects the parameters of a PI-controller such that the Nyquist curve of `P` at the frequency `ωp` is moved to `rl exp(i ϕl)`

If `phasemargin` is supplied, `ϕl` is selected such that the curve is moved to an angle of `phasemargin - 180` degrees

If no `rl` is given, the magnitude of the curve at `ωp` is kept the same and only the phase is affected, the same goes for `ϕl` if no phasemargin is given.

See also `pidplots`, `stabregionPID`
"""
function loopshapingPI(P,ωp; ϕl=0,rl=0, phasemargin = 0, doplot = false)
Pw = P(im*ωp)[1]
ϕp = angle(Pw)
rp = abs(Pw)

if phasemargin > 0
    ϕl = deg2rad(-180+phasemargin)
else
    ϕl = ϕl == 0 ? ϕp : ϕl
end
rl = rl == 0 ? rp : rl

kp = rl/rp*cos(ϕp-ϕl)
ki = rl*ωp/rp*sin(ϕp-ϕl)

C = pid(kp=kp, ki=ki)

if doplot
    gangoffourplot(P,[tf(1),C])
    nyquistplot([P, P*C])
end
return kp,ki,C
end
