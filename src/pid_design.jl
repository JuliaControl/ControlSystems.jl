export pid, pidplots, rlocus, leadlink, laglink, leadlinkat, leadlinkcurve

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
function rlocus(P::LTISystem, K=0)
    K = K == 0 ? linspace(1e-6,500,10000) : K
    Z = tzero(P)
    poles = map(K -> pole(K*P/(1+K*P)), K)
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
    num = [1/a, 1]
    den = [M/a, 1]
    gain = M
    if h <= 0
        return tf(gain*num,den)
    end
    num = c2d_poly2poly(num,h)
    den = c2d_poly2poly(den,h)
    gain *= sum(den)/sum(num)
    return tf(gain*num,den,h)

end


"""
`leadlink(b, N, K; h=0)`

Returns a phase advancing link, the top of the phase curve is located at `ω = b√(N)` where the link amplification is `K√(N)` The bode curve will go from `K`, bend up at `b` and level out at `KN` for frequencies > `bN`

The phase advance at `ω = b√(N)` can be plotted as a function of `N` with `leadlinkcurve()`

See also `leadlinkat`
"""
function leadlink(b, N, K; h=0)
    num = [1/b, 1]
    den = [1/(b*N), 1]
    gain = K
    if h <= 0
        return tf(gain*num,den)
    end
    num = c2d_poly2poly(num,h)
    den = c2d_poly2poly(den,h)
    gain *= sum(den)/sum(num)
    return tf(gain*num,den,h)

end

"""
`leadlinkat(ω, N, K; h=0)`

Returns a phase advancing link, the top of the phase curve is located at `ω` where the link amplification is `K√(N)` The bode curve will go from `K`, bend up at `ω/√(N)` and level out at `KN` for frequencies > `ω√(N)`

The phase advance at `ω` can be plotted as a function of `N` with `leadlinkcurve()`

See also `leadlink`
"""
function leadlinkat(ω, N, K; h=0)
    b = ω / sqrt(N)
    return leadlink(b,N,K,h=h)
end

"""
Plot the phase advance as a function of `N` for a lead link (phase advance link)

See also `Leadlink, leadlinkat`
"""
function leadlinkcurve()
    N = linspace(1,10)
    dph = 180/pi*map(Ni->atan(sqrt(Ni))-atan(1/sqrt(Ni)), N)
    Plots.plot(N,dph, xlabel="N", ylabel="Phase advance [deg]")
end
