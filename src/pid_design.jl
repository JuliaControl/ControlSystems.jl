export pid, pid_tf, pid_ss, pidplots, rlocus, leadlink, laglink, leadlinkat, leadlinkcurve, stabregionPID, loopshapingPI, placePI

"""
    C = pid(param_p, param_i, [param_d]; form=:standard, state_space=false, [Tf], [Ts])

Calculates and returns a PID controller. 

The `form` can be chosen as one of the following
* `:standard` - `Kp*(1 + 1/(Ti*s) + Td*s)` 
* `:series` - `Kc*(1 + 1/(τi*s))*(τd*s + 1)`
* `:parallel` - `Kp + Ki/s + Kd*s`

If `state_space` is set to `true`, either `kd` has to be zero 
or a positive `Tf` has to be provided for creating a filter on 
the input to allow for a state space realization. 
The filter used is `1 / (1 + s*Tf + (s*Tf)^2/2)`, where `Tf` can typically 
be chosen as `Ti/N` for a PI controller and `Td/N` for a PID controller,
and `N` is commonly in the range 2 to 20. 
The state space will be returned on controllable canonical form.

For a discrete controller a positive `Ts` can be supplied.

## Examples
```
C1 = pid(3.3, 1, 2)                             # Kd≠0 works without filter in tf form
C2 = pid(3.3, 1, 2; Tf=0.3, state_space=true)   # In statespace a filter is needed
C3 = pid(2., 3, 0; Ts=0.4, state_space=true)    # Discrete
```

The functions `pid_tf` and `pid_ss` are also exported. They take the same parameters
and is what is actually called in `pid` based on the `state_space` parameter.
"""
function pid(param_p, param_i, param_d=zero(typeof(param_p)); form=:standard, Ts=nothing, Tf=nothing, state_space=false)
    if state_space # Type unstability? Can it be fixed easily, does it matter?
        pid_ss(param_p, param_i, param_d; form, Ts, Tf)
    else
        pid_tf(param_p, param_i, param_d; form, Ts, Tf)
    end
end

function pid_tf(param_p, param_i, param_d=zero(typeof(param_p)); form=:standard, Ts=nothing, Tf=nothing)
    Kp, Ti, Td = convert_pidparams_to_standard(param_p, param_i, param_d, form)
    TE = isnothing(Ts) ? Continuous() : Discrete(Ts)
    if isnothing(Tf)
        if Ti != Inf
            return tf([Kp * Td, Kp, Kp / Ti], [1, 0], TE)
        else
            return tf([Kp * Td, Kp], [1], TE)
        end
    else
        if Ti != Inf
            return tf([Kp * Td, Kp, Kp / Ti], [Tf^2/2, Tf, 1, 0], TE)
        else
            return tf([Kp * Td, Kp], [Tf^2/2, Tf, 1], TE)
        end
    end
end

function pid_ss(param_p, param_i, param_d=zero(typeof(param_p)); form=:standard, Ts=nothing, Tf=nothing)
    Kp, Ti, Td = convert_pidparams_to_standard(param_p, param_i, param_d, form)
    TE = isnothing(Ts) ? Continuous() : Discrete(Ts)
    if !isnothing(Tf)
        A = [0 1 0; 0 0 1; 0 -2/Tf^2 -2/Tf]
        B = [0; 0; 1]
        C = 2 * Kp / Tf^2 * [1/Ti 1 Td]
        D = 0
    elseif Td == 0
        A = 0
        B = 1
        C = Kp / Ti
        D = Kp
    else
        throw(DomainError("cannot create controller as a state space if Td != 0 without a filter. Either create the controller as a transfer function, pid(TransferFunction; params...), or supply Tf to create a filter."))
    end
    return StateSpace(A, B, C, D, TE)
end

"""
    pidplots(P, args...; params_p, params_i, params_d=0, form=:standard, ω=0, grid=false, kwargs...)

Plots interesting figures related to closing the loop around process `P` with a PID controller supplied in `params`
on one of the following forms:
* `:standard` - `Kp*(1 + 1/(Ti*s) + Td*s)` 
* `:series` - `Kc*(1 + 1/(τi*s))*(τd*s + 1)`
* `:parallel` - `Kp + Ki/s + Kd*s`
The sent in values can be arrays to evaluate multiple different controllers, and if `grid=true` it will be a grid search 
over all possible combinations of the values.

Available plots are `:gof` for Gang of four, `:nyquist`, `:controller` for a bode plot of the controller TF and `:pz` for pole-zero maps
and should be supplied as additional arguments to the function.

One can also supply a frequency vector `ω` to be used in Bode and Nyquist plots.

See also `loopshapingPI`, `stabregionPID`
"""
function pidplots(P::LTISystem, args...; 
    params_p, params_i, params_d=0, 
    form=:standard, ω=exp10.(range(-3, stop=3, length=500)), grid=false, 
    kwargs...
)
    if grid
        kps = [i for i in params_p for _ in params_i for _ in params_d]
        kis = [j for _ in params_p for j in params_i for _ in params_d]
        kds = [k for _ in params_p for _ in params_i for k in params_d]
    else
        n = max(length(params_p), length(params_i), length(params_d))
        kps = params_p isa Real ? fill(params_p, n) : params_p
        kis = params_i isa Real ? fill(params_i, n) : params_i
        kds = params_d isa Real ? fill(params_d, n) : params_d
    end

    Cs = LTISystem[]
    PCs = LTISystem[]
    Ts  = LTISystem[]
    labels = Array{String,2}(undef, 1,length(kps))
    for i in eachindex(kps)
        kp = kps[i]
        ki = kis[i]
        kd = kds[i]
        label = latexstring("k_p = $(round(kp, digits=3)),      k_i = $(round(ki, digits=3)),      k_d = $(round(kd, digits=3))")

        C = pid(kp,ki,kd,form=form)
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
    Z = tzeros(P)
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
    laglink(a, M; [Ts])

Returns a phase retarding link, the rule of thumb `a = 0.1ωc` guarantees less than 6 degrees phase margin loss. The bode curve will go from `M`, bend down at `a/M` and level out at 1 for frequencies > `a`
"""
function laglink(a, M; h=nothing, Ts=nothing)
    if !isnothing(h)
        Base.depwarn("`laglink($a, $M; h=$h)` is deprecated, use `laglink($a, $M; Ts=$h)` instead.", Core.Typeof(laglink).name.mt.name)
        Ts = h
    end
    Ts ≥ 0 || throw(ArgumentError("Negative `Ts` is not supported."))
    numerator = [1/a, 1]
    denominator = [M/a, 1]
    gain = M
    G = tf(gain*numerator,denominator)
    return  isnothing(Ts) ? G : c2d(G,Ts)
end


"""
    leadlink(b, N, K; [Ts])

Returns a phase advancing link, the top of the phase curve is located at `ω = b√(N)` where the link amplification is `K√(N)` The bode curve will go from `K`, bend up at `b` and level out at `KN` for frequencies > `bN`

The phase advance at `ω = b√(N)` can be plotted as a function of `N` with `leadlinkcurve()`

Values of `N < 1` will give a phase retarding link.

See also `leadlinkat` `laglink`
"""
function leadlink(b, N, K; h=nothing, Ts=nothing)
    if !isnothing(h)
        Base.depwarn("`leadlink($b, $N, $K; h=$h)` is deprecated, use `leadlink($b, $N, $K; Ts=$h)` instead.", Core.Typeof(leadlink).name.mt.name)
        Ts = h
    end
    Ts ≥ 0 || throw(ArgumentError("Negative `Ts` is not supported."))
    numerator = [1/b, 1]
    denominator = [1/(b*N), 1]
    gain = K
    G = tf(gain*numerator,denominator)
    return  isnothing(Ts) ? G : c2d(G,Ts)

end

"""
    leadlinkat(ω, N, K; [Ts])

Returns a phase advancing link, the top of the phase curve is located at `ω` where the link amplification is `K√(N)` The bode curve will go from `K`, bend up at `ω/√(N)` and level out at `KN` for frequencies > `ω√(N)`

The phase advance at `ω` can be plotted as a function of `N` with `leadlinkcurve()`

Values of `N < 1` will give a phase retarding link.

See also `leadlink` `laglink`
"""
function leadlinkat(ω, N, K; Ts=nothing)
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
    RecipesBase.plot(N,dph, xlabel="N", ylabel="Phase advance [deg]")
end


"""
    fig, kp, ki = stabregionPID(P, [ω]; kd=0)

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
function stabregionPID(P, ω = _default_freq_vector(P,Val{:bode}()); kd=0, form=:standard)
    Pv  = freqresp(P,ω)[:,1,1]
    r   = abs.(Pv)
    phi = angle.(Pv)
    kp  = -cos.(phi)./r
    ki  = kd.*ω.^2 .- ω.*sin.(phi)./r
    plt = RecipesBase.plot(kp,ki,linewidth = 1.5, xlabel=L"k_p", ylabel=L"k_i", title="Stability region of P, k_d = $(round(kd, digits=4))")
    params = convert_pidparams_from_standard(convert_pidparams_to_standard(kp, ki, kd, :parallel)..., form)[1:2]
    plt, params...
end


function stabregionPID(P::Function, ω = exp10.(range(-3, stop=1, length=50)); kd=0, form=:standard)
    Pv      = P.(im*ω)
    r       = abs.(Pv)
    phi     = angle.(Pv)
    kp      = -cos.(phi)./r
    ki      = kd.*ω.^2 .- ω.*sin.(phi)./r
    plt = RecipesBase.plot(kp,ki,linewidth = 1.5, xlabel=L"k_p", ylabel=L"k_i", title="Stability region of P, k_d = $(round(kd, digits=4))")
    params = convert_pidparams_from_standard(convert_pidparams_to_standard(kp, ki, kd, :parallel)..., form)[1:2]
    plt, params...
end


"""
    C, kp, ki = loopshapingPI(P, ωp; ϕl, rl, phasemargin, form=:standard, doplot=false)

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

    kp = rl/rp*cos(ϕp-ϕl)
    ki = rl*ωp/rp*sin(ϕp-ϕl)

    Kp, Ti, Td = convert_pidparams_to_standard(kp, ki, 0, :parallel)
    C = pid(Kp, Ti, Td)

    if doplot
        gangoffourplot(P,[tf(1),C]) |> display
        nyquistplot([P, P*C]) |> display
    end
    C, convert_pidparams_from_standard(Kp, Ti, Td, form)...
end


"""
    C, params = placePI(P, ω₀, ζ; form=:standard)

Selects the parameters of a PI-controller such that the poles of 
closed loop between `P` and `C` are placed to match the poles of 
`s^2 + 2ζω₀s + ω₀^2`.

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
    kp = -tmp / (a^2*ω₀^2 - 2*a*b*ω₀*ζ + b^2)
    ki = tmp / (ω₀^2*(a*d - b*c))
    pid(kp, ki), convert_pidparams_from_standard(kp, ki, 0, form)[1:2]...
end

placePI(sys::LTISystem, args...; kwargs...) = placePI(tf(sys), args...; kwargs...)

"""
    Kp, Ti, Td = convert_pidparams_to_standard(param_p, param_i, param_d, form)

Convert parameters from form `form` to `:standard` form. 

The `form` can be chosen as one of the following
* `:standard` - `Kp*(1 + 1/(Ti*s) + Td*s)` 
* `:series` - `Kc*(1 + 1/(τi*s))*(τd*s + 1)`
* `:parallel` - `Kp + Ki/s + Kd*s`
"""
function convert_pidparams_to_standard(param_p, param_i, param_d, form)
    if form === :standard
        return param_p, param_i, param_d
    elseif form === :series
        return (
            param_p * (param_i + param_d) / param_i,
            param_i + param_d,
            param_i * param_d / (param_i + param_d)
        )
    elseif form === :parallel
        return (param_p, param_p / param_i, param_d / param_p)
    else
        throw(ArgumentError("form $(form) not supported."))
    end
end

"""
    param_p, param_i, param_d = convert_pidparams_from_standard(Kp, Ti, Td, form)

Convert parameters to form `form` from `:standard` form. 

The `form` can be chosen as one of the following
* `:standard` - `Kp*(1 + 1/(Ti*s) + Td*s)` 
* `:series` - `Kc*(1 + 1/(τi*s))*(τd*s + 1)`
* `:parallel` - `Kp + Ki/s + Kd*s`
"""
function convert_pidparams_from_standard(Kp, Ti, Td, form)
    if form === :standard
        return Kp, Ti, Td
    elseif form === :series
        return (
            (Ti - sqrt(Ti * (Ti - 4 * Td))) / 2 * Kp / Ti,
            (Ti - sqrt(Ti * (Ti - 4 * Td))) / 2,
            (Ti + sqrt(Ti * (Ti - 4 * Td))) / 2,
        )
    elseif form === :parallel
        return (Kp, Kp/Ti, Td*Kp)
    else
        throw(ArgumentError("form $(form) not supported."))
    end
end
