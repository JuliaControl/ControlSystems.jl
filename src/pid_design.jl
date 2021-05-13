export pid, pidplots, rlocus, leadlink, laglink, leadlinkat, leadlinkcurve, stabregionPID, loopshapingPI

"""
    C = pid(; kp=0, ki=0; kd=0, time=false, series=false)

Calculates and returns a PID controller on transfer function form.
- `time` indicates whether or not the parameters are given as gains (default) or as time constants
- `series` indicates  whether or not the series form or parallel form (default) is desired
"""
function pid(; kp=0., ki=0., kd=0., time=false, series=false)
    s = tf("s")
    if series
        return time ? kp*(one(kp) + one(kp)/(ki*s) + kd*s) : kp*(one(kp) + ki/s + kd*s)
    else
        return time ? kp + one(kp)/(ki*s) + kd*s : kp + ki/s + kd*s
    end
end

"""
    pidplots(P, args...; kps=0, kis=0, kds=0, time=false, series=false, ω=0)

Plots interesting figures related to closing the loop around process `P` with a PID controller
Send in a bunch of PID-parameters in any of the vectors kp, ki, kd. The vectors must be the same length.

-`time` indicates whether or not the parameters are given as gains (default) or as time constants
-`series` indicates  whether or not the series form or parallel form (default) is desired

Available plots are `:gof` for Gang of four, `:nyquist`, `:controller` for a bode plot of the controller TF and `:pz` for pole-zero maps

One can also supply a frequency vector ω to be used in Bode and Nyquist plots

See also `loopshapingPI`, `stabregionPID`
"""
function pidplots(P::LTISystem, args...; kps=0, kis=0, kds=0, time=false, series=false, ω=0, grid = false, kwargs...)
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
    ω   = ω   == 0 ? exp10.(range(-3, stop=3, length=500)) : ω

    getColorSys(i)   = convert(Colors.RGB,Colors.HSV(360*((i-1)/(length(kps)))^1.5,0.9,0.8))

    Cs = LTISystem[]
    PCs = LTISystem[]
    Ts  = LTISystem[]
    labels = Array{String,2}(undef, 1,length(kps))
    colors =  Array{Colors.RGB{Float64},2}(undef, 1, length(kps))
    for (i,kp) = enumerate(kps)
        ki = kis[i]
        kd = kds[i]
        label = latexstring("k_p = $(round(kp, digits=3)),      k_i = $(round(ki, digits=3)),      k_d = $(round(kd, digits=3))")

        C = pid(kp=kp,ki=ki,kd=kd,time=time,series=series)
        S,D,N,T = gangoffour(P,C)
        push!(Cs, C)
        push!(PCs, P*C)
        push!(Ts, T)
        labels[i] = label
    end

    if :nyquist ∈ args
        nyquistplot(PCs, ω; lab=labels, title="Nyquist curves", kwargs...) |> display
    end
    if :gof ∈ args
        gangoffourplot(P, Cs, ω; lab=labels, kwargs...) |> display
    end
    if :pz ∈ args
        pzmap(Ts; title="Pole-zero map", kwargs...) |> display
    end
    if :controller ∈ args
        bodeplot(Cs, ω; lab=labels, title="Controller bode plot", kwargs...) |> display
    end
end

@userplot Rlocusplot
@deprecate rlocus(args...;kwargs...) rlocusplot(args...;kwargs...)



function getpoles(G, K) # If OrdinaryDiffEq is installed, we override getpoles with an adaptive method
    P          = numpoly(G)[1]
    Q          = denpoly(G)[1]
    f          = (y,_,k) -> sort(ComplexF64.(Polynomials.roots(k[1]*P+Q)), by=imag)
    prob       = OrdinaryDiffEq.ODEProblem(f,f(0.,0.,0.),(0.,K[end]))
    integrator = OrdinaryDiffEq.init(prob,OrdinaryDiffEq.Tsit5(),reltol=1e-8,abstol=1e-8)
    ts         = Vector{Float64}()
    poleout    = Vector{Vector{ComplexF64}}()
    for i in integrator
       push!(poleout,integrator.k[1])
       push!(ts,integrator.t[1])
    end
    poleout = hcat(poleout...)'
    poleout, ts
end



"""
    rlocusplot(P::LTISystem, K)

Computes and plots the root locus of the SISO LTISystem P with
a negative feedback loop and feedback gains `K`, if `K` is not provided,
range(1e-6,stop=500,length=10000) is used.
If `OrdinaryDiffEq.jl` is installed and loaded by the user (`using OrdinaryDiffEq`), `rlocusplot` will use an adaptive step-size algorithm to
select values of `K`. A scalar `Kmax` can then be given as second argument.
"""
rlocus
@recipe function rlocus(p::Rlocusplot; K=500)
    P = p.args[1]
    K = K isa Number ? range(1e-6,stop=K,length=10000) : K
    Z = tzero(P)
    poles, K = getpoles(P,K)
    redata = real.(poles)
    imdata = imag.(poles)
    ylim = (max(-50,minimum(imdata)), min(50,maximum(imdata)))
    xlim = (max(-50,minimum(redata)), min(50,maximum(redata)))
    framestyle --> :zerolines
    title --> "Root locus"
    xguide --> "Re(roots)"
    yguide --> "Im(roots)"
    form(k, p) = Printf.@sprintf("%.4f", k) * "  pole=" * Printf.@sprintf("%.3f%+.3fim", real(p), imag(p))
    @series begin
        legend --> false
        ylims  --> ylim
        xlims  --> xlim
        hover := "K=" .* form.(K,poles)
        label := ""
        redata, imdata
    end
    @series begin
        seriestype := :scatter
        markershape --> :circle
        markersize --> 10
        label --> "Zeros"
        real.(Z), imag.(Z)
    end
    @series begin
        seriestype := :scatter
        markershape --> :xcross
        markersize --> 10
        label --> "Open-loop poles"
        redata[1,:], imdata[1,:]
    end
end

"""
    laglink(a, M; Ts=0)

Returns a phase retarding link, the rule of thumb `a = 0.1ωc` guarantees less than 6 degrees phase margin loss. The bode curve will go from `M`, bend down at `a/M` and level out at 1 for frequencies > `a`
"""
function laglink(a, M; h=nothing, Ts=0)
    if !isnothing(h)
        Base.depwarn("`laglink($a, $M; h=$h)` is deprecated, use `laglink($a, $M; Ts=$h)` instead.", Core.Typeof(laglink).name.mt.name)
        Ts = h
    end
    @assert Ts ≥ 0 "Negative `Ts` is not supported."
    numerator = [1/a, 1]
    denominator = [M/a, 1]
    gain = M
    G = tf(gain*numerator,denominator)
    return  Ts <= 0 ? G : c2d(G,Ts)
end


"""
    leadlink(b, N, K; Ts=0)

Returns a phase advancing link, the top of the phase curve is located at `ω = b√(N)` where the link amplification is `K√(N)` The bode curve will go from `K`, bend up at `b` and level out at `KN` for frequencies > `bN`

The phase advance at `ω = b√(N)` can be plotted as a function of `N` with `leadlinkcurve()`

Values of `N < 1` will give a phase retarding link.

See also `leadlinkat` `laglink`
"""
function leadlink(b, N, K; h=nothing, Ts=0)
    if !isnothing(h)
        Base.depwarn("`leadlink($b, $N, $K; h=$h)` is deprecated, use `leadlink($b, $N, $K; Ts=$h)` instead.", Core.Typeof(leadlink).name.mt.name)
        Ts = h
    end
    @assert Ts ≥ 0 "Negative `Ts` is not supported."
    numerator = [1/b, 1]
    denominator = [1/(b*N), 1]
    gain = K
    G = tf(gain*numerator,denominator)
    return  Ts <= 0 ? G : c2d(G,Ts)

end

"""
    leadlinkat(ω, N, K; Ts=0)

Returns a phase advancing link, the top of the phase curve is located at `ω` where the link amplification is `K√(N)` The bode curve will go from `K`, bend up at `ω/√(N)` and level out at `KN` for frequencies > `ω√(N)`

The phase advance at `ω` can be plotted as a function of `N` with `leadlinkcurve()`

Values of `N < 1` will give a phase retarding link.

See also `leadlink` `laglink`
"""
function leadlinkat(ω, N, K; h=nothing, Ts=0)
    if !isnothing(h)
        Base.depwarn("`leadlinkat($ω, $N, $K; h=$h)` is deprecated, use `leadlinkat($ω, $N, $K; Ts=$h)` instead.", Core.Typeof(leadlinkat).name.mt.name)
        Ts = h
    end
    b = ω / sqrt(N)
    return leadlink(b,N,K,Ts=Ts)
end

"""
leadlinkcurve(start=1)

Plot the phase advance as a function of `N` for a lead link (phase advance link)
If an input argument `s` is given, the curve is plotted from `s` to 10, else from 1 to 10.

See also `Leadlink, leadlinkat`
"""
function leadlinkcurve(start=1)
    N = range(start, stop=10, length=50)
    dph = 180/pi*map(Ni->atan(sqrt(Ni))-atan(1/sqrt(Ni)), N)
    Plots.plot(N,dph, xlabel="N", ylabel="Phase advance [deg]")
end



"""
    fig, kp, ki = stabregionPID(P, [ω]; kd=0, doplot = true)

Segments of the curve generated by this program
is the boundary of the stability region for a
process with transfer function P(s)
The PID controller is assumed to be on the form kp +ki/s +kd s

The curve is found by analyzing
P(s)*C(s) = -1 ⟹
|PC| = |P| |C| = 1
arg(P) + arg(C) = -π

If `P` is a function (e.g. s -> exp(-sqrt(s)) ), the stability of feedback loops using PI-controllers can be analyzed for processes with models with arbitrary analytic functions

See also `stabregionPID`, `loopshapingPI`, `pidplots`
"""
function stabregionPID(P, ω = _default_freq_vector(P,Val{:bode}()); kd=0, doplot = true)
    Pv  = freqresp(P,ω)[:,1,1]
    r   = abs.(Pv)
    phi = angle.(Pv)
    kp  = -cos.(phi)./r
    ki  = kd.*ω.^2 .- ω.*sin.(phi)./r
    Plots.plot(kp,ki,linewidth = 1.5, xlabel=L"k_p", ylabel=L"k_i", title="Stability region of P, k_d = $(round(kd, digits=4))"), kp, ki
end


function stabregionPID(P::Function, ω = exp10.(range(-3, stop=1, length=50)); kd=0, doplot = true)
    Pv      = P.(im*ω)
    r       = abs.(Pv)
    phi     = angle.(Pv)
    kp      = -cos.(phi)./r
    ki      = kd.*ω.^2 .- ω.*sin.(phi)./r
    Plots.plot(kp,ki,linewidth = 1.5, xlabel=L"k_p", ylabel=L"k_i", title="Stability region of P, k_d = $(round(kd, digits=4))"), kp, ki
end


"""
    kp,ki,C = loopshapingPI(P,ωp; ϕl,rl, phasemargin, doplot = false)

Selects the parameters of a PI-controller such that the Nyquist curve of `P` at the frequency `ωp` is moved to `rl exp(i ϕl)`

If `phasemargin` is supplied, `ϕl` is selected such that the curve is moved to an angle of `phasemargin - 180` degrees

If no `rl` is given, the magnitude of the curve at `ωp` is kept the same and only the phase is affected, the same goes for `ϕl` if no phasemargin is given.

Set `doplot = true` to plot the `gangoffourplot` and `nyquistplot` of the system.

See also `pidplots`, `stabregionPID`
"""
function loopshapingPI(P,ωp; ϕl=0,rl=0, phasemargin = 0, doplot = false)
    Pw = P(im*ωp)[1]
    ϕp = angle(Pw)
    rp = abs.(Pw)

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
        gangoffourplot(P,[tf(1),C]) |> display
        nyquistplot([P, P*C]) |> display
    end
    return kp,ki,C
end
