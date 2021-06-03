export pid, pidplots, rlocus, leadlink, laglink, leadlinkat, leadlinkcurve, stabregionPID, loopshapingPI, placePI

"""
    C = pid([type=StateSpace;] Tf=0, params...)

Calculates and returns a PID controller with `type` representation.
Default is `StateSpace`, but the first argument can be supplied as
`TransferFunction` to create it as a transfer function.

If state space is selected, either the derivative part has to be zero 
or `Tf` has to be provided for creating a filter on the input to 
allow for a state space realization. 
The filter used is `1 / (1 + s*Tf + (s*Tf)^2/2)`, where `Tf` is normally
chosen as `Ti/N` for a PI controller and `Td/N` for a PID controller,
and `N` is commonly in the range 2 to 20.

The state space will be returned on controllable canonical form.

`params` are supplied as keyword arguments with names corresponding
to either one of the following forms:
* `:standard` - `Kp*(1 + 1/(Ti*s) + Td*s)` 
* `:series` - `Kc*(1 + 1/(τi*s))*(1 + τd*s)`
* `:parallel` - `Kp + Ki/s + Kd*s`
Default values if only a subset of a form is supplied is `Kp=Kc=1`, 
`Ti=τi=Inf`, `Ki=Kd=Td=τd=0`.

## Examples
```
C1 = pid(form=:parallel, Kd=4, Ki=2, Kd=1)
C2 = pid(Kp=2, Ti=3, Td=1, N=4)
C3 = pid(StateSpace; Kp=2, Ti=3)
C4 = pid(StateSpace; form=:standard, Kp=2, Ti=3, Td=1, N=4)
```
"""
pid(; kwargs...) = pid(StateSpace; kwargs...)

function pid(::Type{TransferFunction}; Tf=0, params...)
    p = convert_pid_params(:standard; params...)
    if p.Ti != Inf
        tf([p.Kp * p.Td, p.Kp, p.Kp / p.Ti], [Tf^2/2, Tf, 1, 0])
    else
        tf([p.Kp * p.Td, p.Kp], [Tf^2/2, Tf, 1])
    end
end

function pid(::Type{StateSpace}; Tf=0, params...)
    p = convert_pid_params(:standard; params...)
    if Tf != 0
        A = [0 1 0; 0 0 1; 0 -2/Tf^2 -2/Tf]
        B = [0; 0; 1]
        C = 2 * p.Kp / Tf^2 * [1/p.Ti 1 p.Td]
        D = 0
    elseif p.Td == 0
        A = 0
        B = 1
        C = p.Kp / p.Ti
        D = p.Kp
    else
        throw(DomainError("cannot create a state space form for Td != 0 and Tf == 0."))
    end
    return ss(A, B, C, D)
end

"""
    pidplots(P, args...; ω=0, grid=false, params..., kwargs...)

Plots interesting figures related to closing the loop around process `P` with a PID controller supplied in `params`
on one of the following forms:
* `:standard` - `Kp*(1 + 1/(Ti*s) + Td*s)` 
* `:series` - `Kc*(1 + 1/(τi*s))*(1 + τd*s)`
* `:parallel` - `Kp + Ki/s + Kd*s`
The sent in values can be arrays to evaluate multiple different controllers, and if `grid=true` it will be a grid search 
over all possible combinations of the values.

Available plots are `:gof` for Gang of four, `:nyquist`, `:controller` for a bode plot of the controller TF and `:pz` for pole-zero maps
and should be supplied as additional arguments to the function.

One can also supply a frequency vector `ω` to be used in Bode and Nyquist plots.

See also `loopshapingPI`, `stabregionPID`
"""
function pidplots(P::LTISystem, args...; 
        Kp=0.0, Kc=0.0, Ki=0.0, Kd=0.0, Ti=eltype(Kp)(Inf), τi=eltype(Kc)(Inf), Td=0.0, τd=0.0, 
        ω=0, grid = false, kwargs...)
    params = if Kc != 0
        (;Kc, τi, τd)
    elseif Ki != 0 || Kd != 0
        (;Kp, Ki, Kd)
    elseif Kp != 0
        (;Kp, Ti, Td)
    else
        throw(ArgumentError("the supplied parameters does not sufficiently describe a supported form."))
    end
    kps, kis, kds = values(params)

    if grid
        kp = [i for i in kps, _ in kis, _ in kds][:]
        ki = [j for _ in kps, j in kis, _ in kds][:]
        kd = [k for _ in kps, _ in kis, k in kds][:]
        kps, kis, kds = kp, ki, kd
    else
        n = max(length(kps), length(kis), length(kds))
        kps = kps isa Number ? fill(kps, n) : kps
        kis = kis isa Number ? fill(kis, n) : kis
        kds = kds isa Number ? fill(kds, n) : kds
    end
    ω = ω == 0 ? exp10.(range(-3, stop=3, length=500)) : ω

    getColorSys(i) = convert(Colors.RGB,Colors.HSV(360*((i-1)/(length(kps)))^1.5,0.9,0.8))

    Cs = LTISystem[]
    PCs = LTISystem[]
    Ts  = LTISystem[]
    labels = Array{String,2}(undef, 1,length(kps))
    for i = eachindex(kps)
        param = zip(keys(params), (kps[i], kis[i], kds[i]))
        label = join(("$(k) = $(v)" for (k, v) in param), ", ")

        C = pid(TransferFunction; param...)
        T = robust_minreal(feedback(P*C, 1))
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
    fig, params = stabregionPID(P, [ω]; kd=0, form=:standard)

Segments of the curve generated by this program
is the boundary of the stability region for a
process with transfer function P(s)

The PI controller is returned on one of the following forms
selected by `form`
* `:standard` - `Kp*(1 + 1/(Ti*s))` 
* `:series` - `Kc*(1 + 1/(τi*s))`
* `:parallel` - `Kp + Ki/s`

The curve is found by analyzing
P(s)*C(s) = -1 ⟹
|PC| = |P| |C| = 1
arg(P) + arg(C) = -π

If `P` is a function (e.g. s -> exp(-sqrt(s)) ), the stability of feedback loops using PI-controllers can be analyzed for processes with models with arbitrary analytic functions

See also `loopshapingPI`, `pidplots`
"""
function stabregionPID(P, ω = _default_freq_vector(P,Val{:bode}()); form=:standard, kd=0)
    Pv  = freqresp(P,ω)[:,1,1]
    r   = abs.(Pv)
    phi = angle.(Pv)
    kp  = -cos.(phi)./r
    ki  = kd.*ω.^2 .- ω.*sin.(phi)./r
    params = [convert_pid_params(form; Kp=kkp, Ki=kki) for (kkp, kki) in zip(kp, ki)]
    kps = map(x->values(x)[1], params)
    kis = map(x->values(x)[2], params)
    params = (; zip(keys(first(params))[1:2], [kps, kis])...)
    Plots.plot(kp,ki,linewidth = 1.5, xlabel=L"k_p", ylabel=L"k_i", title="Stability region of P, k_d = $(round(kd, digits=4))"), params
end


function stabregionPID(P::Function, ω = exp10.(range(-3, stop=1, length=50)); form=:standard, kd=0)
    Pv      = P.(im*ω)
    r       = abs.(Pv)
    phi     = angle.(Pv)
    kp      = -cos.(phi)./r
    ki      = kd.*ω.^2 .- ω.*sin.(phi)./r
    params = [convert_pid_params(form; Kp=kkp, Ki=kki) for (kkp, kki) in zip(kp, ki)]
    kps = map(x->values(x)[1], params)
    kis = map(x->values(x)[2], params)
    params = (; zip(keys(first(params))[1:2], [kps, kis])...)
    Plots.plot(kp,ki,linewidth = 1.5, xlabel=L"k_p", ylabel=L"k_i", title="Stability region of P, k_d = $(round(kd, digits=4))"), params
end


"""
    C, params = loopshapingPI(P, ωp; ϕl, rl, phasemargin, form=:standard, doplot=false)

Selects the parameters of a PI-controller such that the Nyquist curve of `P` at the frequency `ωp` is moved to `rl exp(i ϕl)`

The parameters can be returned as one of several common representations 
chosen by `form`, the options are
* `:standard` - `Kp*(1 + 1/(Ti*s) + Td*s)` 
* `:series` - `Kc*(1 + 1/(τi*s))*(τd*s + 1)`
* `:parallel` - `Kp + Ki/s + Kd*s`

If `phasemargin` is supplied, `ϕl` is selected such that the curve is moved to an angle of `phasemargin - 180` degrees

If no `rl` is given, the magnitude of the curve at `ωp` is kept the same and only the phase is affected, the same goes for `ϕl` if no phasemargin is given.

Set `doplot = true` to plot the `gangoffourplot` and `nyquistplot` of the system.

See also `pidplots`, `stabregionPID`
"""
function loopshapingPI(P, ωp; ϕl=0, rl=0, phasemargin=0, form=:standard, doplot=false)
    Pw = P(im*ωp)[1]
    ϕp = angle(Pw)
    rp = abs.(Pw)

    if phasemargin > 0
        ϕl = deg2rad(-180+phasemargin)
    else
        ϕl = ϕl == 0 ? ϕp : ϕl
    end
    rl = rl == 0 ? rp : rl

    Kp = rl/rp*cos(ϕp-ϕl)
    Ki = rl*ωp/rp*sin(ϕp-ϕl)

    C = pid(;Kp, Ki)

    if doplot
        gangoffourplot(P,[tf(1),C]) |> display
        nyquistplot([P, P*C]) |> display
    end
    return C, convert_pid_params(form; Kp, Ki)
end

"""
    C, params = placePI(P, ω₀, ζ; form=:standard)

Selects the parameters of a PI-controller such that the poles of 
closed loop between `P` and `C` are placed to match the poles of 
`s^2 + 2ζω₀ + ω₀^2`.

The parameters can be returned as one of several common representations 
chose by `form`, the options are
* `:standard` - `Kp*(1 + 1/(Ti*s) + Td*s)` 
* `:series` - `Kc*(1 + 1/(τi*s))*(τd*s + 1)`
* `:parallel` - `Kp + Ki/s + Kd*s`

`C` is the returned transfer function of the controller and `params` 
is a named tuple containing the parameters.
The parameters can be accessed as `params.Kp` or `params["Kp"]` from the named tuple,
or they can be unpacked using `Kp, Ti, Td = values(params)`.
"""
function placePI(P::TransferFunction{<:Continuous, <:SisoRational{T}}, ω₀, ζ; form=:standard) where T
    num = numvec(P)[]
    den = denvec(P)[]
    if length(den) != 2 || length(num) > 2
        throw(DomainError("can only place poles using PI for proper first-order systems"))
    end
    if length(num) == 1
        num = [0; num]
    end
    a, b = num
    c, d = den
    # Calculates PI on standard/series form
    tmp = (a*c*ω₀^2 - 2*b*c*ζ*ω₀ + b*d)
    Kp = -tmp / (a^2*ω₀^2 - 2*a*b*ω₀*ζ + b^2)
    Ti = tmp / (ω₀^2*(a*d - b*c))
    return pid(;Kp, Ti), convert_pid_params(form; Kp, Ti)
end

placePI(sys::LTISystem, args...; kwargs...) = placePI(tf(sys), args...; kwargs...)

"""
    params = convert_pidparams(target; params...)

Convert parameters from a form idendified by the parameter names sent in to `target` form.

`target` can be chosen as one of the following forms
* `:standard` - `Kp*(1 + 1/(Ti*s) + Td*s)` 
* `:series` - `Kc*(1 + 1/(τi*s))*(τd*s + 1)`
* `:parallel` - `Kp + Ki/s + Kd*s`
with default `params` values `Kp=Kc=1`, `Ti=τi=Inf`, `Ki=Kd=Td=τd=0`.

`params` should be supplied with parameter names corresponding to the names used in the above 
equations.
"""
function convert_pid_params(target; params...)
    if haskey(params, :Kc)
        Kc = params[:Kc]
        τi = get(params, :τi, typeof(Kc)(Inf))
        τd = get(params, :τd, 0.0)
        Kp = Kc * (1 + τd / τi)
        Ti = τd + τi
        Td = τd * τi / (τd + τi)
    elseif haskey(params, :Ki) || haskey(params, :Kd)
        Kp = get(params, :Kp, 0.0)
        Ki = get(params, :Ki, 0.0)
        Kd = get(params, :Kd, 0.0)
        Ti = Ki == 0 ? typeof(Kp)(Inf) : Kp / Ki
        Td = Kd == 0 ? 0 : Kd / Kp
    elseif haskey(params, :Kp)
        Kp = get(params, :Kp, 1.0)
        Ti = get(params, :Ti, typeof(Kp)(Inf))
        Td = get(params, :Td, 0.0)
    else
        throw(ArgumentError("the supplied parameters does not sufficiently describe a supported form, params=$(params)"))
    end

    if target === :series
        1 < 4*Td/Ti && throw(DomainError("series form cannot be used with complex zeros."))
        return (
            Kc = Kp/2 * (1 + sqrt(1 - 4*Td/Ti)), 
            τi = Ti/2 * (1 + sqrt(1 - 4*Td/Ti)), 
            τd = Ti == Inf ? Td : Ti/2 * (1 - sqrt(1 - 4*Td/Ti))
        )
    elseif target === :parallel
        return (Kp=Kp, Ki=Kp/Ti, Kd=Kp*Td)
    elseif target === :standard
        return (;Kp, Ti, Td)
    else
        throw(ArgumentError("form $(target) not supported"))
    end
end
