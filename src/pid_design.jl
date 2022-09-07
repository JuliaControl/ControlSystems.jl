export pid, pid_tf, pid_ss, pidplots, rlocus, leadlink, laglink, leadlinkat, leadlinkcurve, stabregionPID, loopshapingPI, placePI, loopshapingPID

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

@deprecate pid(; kp=0, ki=0, kd=0, series = false) pid(kp, ki, kd; form=series ? :series : :parallel)

function pid_tf(param_p, param_i, param_d=zero(typeof(param_p)); form=:standard, Ts=nothing, Tf=nothing)
    Kp, Ti, Td = convert_pidparams_to_standard(param_p, param_i, param_d, form)
    TE = isnothing(Ts) ? Continuous() : Discrete(Ts)
    ia = Ti != Inf && Ti != 0 # integral action, 0 would result in division by zero, but typically indicates that the user wants no integral action
    if isnothing(Tf)
        if ia
            return tf([Kp * Td, Kp, Kp / Ti], [1, 0], TE)
        else
            return tf([Kp * Td, Kp], [1], TE)
        end
    else
        if ia
            return tf([Kp * Td, Kp, Kp / Ti], [Tf^2/2, Tf, 1, 0], TE)
        else
            return tf([Kp * Td, Kp], [Tf^2/2, Tf, 1], TE)
        end
    end
end

function pid_ss(param_p, param_i, param_d=zero(typeof(param_p)); form=:standard, Ts=nothing, Tf=nothing)
    Kp, Ti, Td = convert_pidparams_to_standard(param_p, param_i, param_d, form)
    TE = isnothing(Ts) ? Continuous() : Discrete(Ts)
    ia = Ti != Inf && Ti != 0 # integral action, 0 would result in division by zero, but typically indicates that the user wants no integral action
    if !isnothing(Tf)
        if ia
            A = [0 1 0; 0 0 1; 0 -2/Tf^2 -2/Tf]
            B = [0; 0; 1]
            C = 2 * Kp / Tf^2 * [1/Ti 1 Td]
        else
            A = [0 1; -2/Tf^2 -2/Tf]
            B = [0; 1]
            C = 2 * Kp / Tf^2 * [1 Td]
        end
        D = 0
    elseif Td == 0
        if ia
            A = 0
            B = 1
            C = Kp / Ti # Ti == 0 would result in division by zero, but typically indicates that the user wants no integral action
            D = Kp
        else
            return StateSpace([Kp], TE)
        end
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


function getpoles(G, K)
    issiso(G) || error("root locus only supports SISO systems")
    G isa TransferFunction || (G = tf(G))
    P = numpoly(G)[]
    Q = denpoly(G)[]
    T = float(eltype(K))
    ϵ = eps(T)
    nx = length(Q)
    D = zeros(nx-1, nx-1) # distance matrix
    prevpoles = ComplexF64[]
    temppoles = zeros(ComplexF64, nx-1)
    f = function (y,_,k)
        if k == 0 && length(P) > length(Q) 
            # More zeros than poles, make sure the vector of roots is of correct length when k = 0
            # When this happens, there are fewer poles for k = 0, these poles can be seen as beeing located somewhere at Inf
            # We get around the problem by not allowing k = 0 for non-proper systems.
            k = ϵ
        end
        newpoles = ComplexF64.(Polynomials.roots(k[1]*P+Q))
        if !isempty(prevpoles)
            D .= abs.(newpoles .- transpose(prevpoles))
            assignment, cost = Hungarian.hungarian(D)
            for i = 1:nx-1
                temppoles[assignment[i]] = newpoles[i]
            end
            newpoles .= temppoles
        end
        prevpoles = newpoles
        newpoles
    end
    prob       = OrdinaryDiffEq.ODEProblem(f,f(0.,0.,0),(0,K[end]))
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
    rlocusplot(P::LTISystem; K)

Computes and plots the root locus of the SISO LTISystem P with
a negative feedback loop and feedback gains between 0 and `K`. `rlocusplot` will use an adaptive step-size algorithm to
determine the values of the feedback gains used to generate the plot.
"""
rlocusplot
@recipe function rlocusplot(p::Rlocusplot; K=500)
    P = p.args[1]
    K = K isa Number ? range(1e-6,stop=K,length=10000) : K
    Z = tzeros(P)
    roots, K = getpoles(P,K)
    redata = real.(roots)
    imdata = imag.(roots)
    
    ylims --> (max(-50,minimum(imdata) - 1), min(50,maximum(imdata) + 1))
    xlims --> (max(-50,minimum(redata) - 1), clamp(maximum(redata) + 1, 1, 50))
    framestyle --> :zerolines
    title --> "Root locus"
    xguide --> "Re(roots)"
    yguide --> "Im(roots)"
    form(k, p) = Printf.@sprintf("%.4f", k) * "  pole=" * Printf.@sprintf("%.3f%+.3fim", real(p), imag(p))
    @series begin
        legend --> false
        hover := "K=" .* form.(K,roots)
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

```math
\\dfrac{s + a}{s + a/M} = M \\dfrac{1 + s/a}{1 + sM/a}
```
"""
function laglink(a, M; Ts=nothing)
    Ts !== nothing && (Ts ≥ 0 || throw(ArgumentError("Negative `Ts` is not supported.")))
    M > 1 || @warn "M should be ≥ 1 for the link to be phase retarding (increasing gain)"
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

```math
KN \\dfrac{s + b}{s + bN} = K \\dfrac{1 + s/b}{1 + s/(bN)}
```

See also `leadlinkat` `laglink`
"""
function leadlink(b, N, K; h=nothing, Ts=nothing)
    Ts !== nothing && (Ts ≥ 0 || throw(ArgumentError("Negative `Ts` is not supported.")))
    N > 1 || @warn "N should be ≥ 1 for the link to be phase advancing."
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

@userplot LeadLinkCurve

"""
    leadlinkcurve(start=1)

Plot the phase advance as a function of `N` for a lead link (phase advance link)
If an input argument `start` is given, the curve is plotted from `start` to 10, else from 1 to 10.

See also `leadlink, leadlinkat`
"""
leadlinkcurve
@recipe function leadlinkcurve(p::LeadLinkCurve)
    start = isempty(p.args) ? 1 : p.args[1]
    N = range(start, stop=10, length=50)
    dph = map(Ni->180/pi*atan(sqrt(Ni))-atan(1/sqrt(Ni)), N)
    @series begin
        xlabel := "N"
        ylabel := "Phase advance [deg]"
        N,dph
    end
end


"""
    kp, ki, fig = stabregionPID(P, [ω]; kd=0, doplot=false, form=:standard)

Segments of the curve generated by this program
is the boundary of the stability region for a
process with transfer function P(s)
The provided derivative gain is expected on parallel form, i.e., the form kp + ki/s + kd s, but the result can be transformed to any form given by the `form` keyword.
The curve is found by analyzing
```math
P(s)C(s) = -1 ⟹ \\\\
|PC| = |P| |C| = 1 \\\\
arg(P) + arg(C) = -π
```
If `P` is a function (e.g. s -> exp(-sqrt(s)) ), the stability of feedback loops using PI-controllers can be analyzed for processes with models with arbitrary analytic functions
See also [`loopshapingPI`](@ref), [`loopshapingPID`](@ref), [`pidplots`](@ref)
"""
function stabregionPID(P, ω = _default_freq_vector(P,Val{:bode}()); kd=0, form=:standard, doplot=false)
    Pv  = freqrespv(P,ω)
    r   = abs.(Pv)
    phi = angle.(Pv)
    kp  = @. -cos(phi)/r
    ki  = @. kd*ω^2 - ω*sin(phi)/r
    kp, ki  = convert_pidparams_from_to(kp, ki, kd, :parallel, form)
    fig = if doplot
        RecipesBase.plot(kp,ki,linewidth = 1.5, xlabel=L"k_p", ylabel=L"k_i", title="Stability region of P, k_d = $(round(kd, digits=4))")
    else 
        nothing
    end
    kp, ki, fig
end


function stabregionPID(P::Function, ω = exp10.(range(-3, stop=1, length=50)); kd=0, form=:standard, doplot=false)
    Pv      = P.(im*ω)
    r       = abs.(Pv)
    phi     = angle.(Pv)
    kp      = -cos.(phi)./r
    ki      = @. kd*ω^2 - ω*sin(phi)/r
    kp, ki  = convert_pidparams_from_to(kp, ki, kd, :parallel, form)
    fig = if doplot
        RecipesBase.plot(kp,ki,linewidth = 1.5, xlabel=L"k_p", ylabel=L"k_i", title="Stability region of P, k_d = $(round(kd, digits=4))")
    else 
        nothing
    end
    kp, ki, fig
end


"""
    C, kp, ki, fig, CF = loopshapingPI(P, ωp; ϕl, rl, phasemargin, form=:standard, doplot=false, Tf, F)

Selects the parameters of a PI-controller (on parallel form) such that the Nyquist curve of `P` at the frequency `ωp` is moved to `rl exp(i ϕl)`

The parameters can be returned as one of several common representations 
chosen by `form`, the options are
* `:standard` - ``K_p(1 + 1/(T_i s) + T_ds)``
* `:series` - ``K_c(1 + 1/(τ_i s))(τ_d s + 1)``
* `:parallel` - ``K_p + K_i/s + K_d s``

If `phasemargin` is supplied (in degrees), `ϕl` is selected such that the curve is moved to an angle of `phasemargin - 180` degrees

If no `rl` is given, the magnitude of the curve at `ωp` is kept the same and only the phase is affected, the same goes for `ϕl` if no phasemargin is given.

- `Tf`: An optional time constant for second-order measurement noise filter on the form `tf(1, [Tf^2, 2*Tf/sqrt(2), 1])` to make the controller strictly proper.
- `F`: A pre-designed filter to use instead of the default second-order filter that is used if `Tf` is given.
- `doplot` plot the `gangoffourplot` and `nyquistplot` of the system.

See also [`loopshapingPID`](@ref), [`pidplots`](@ref), [`stabregionPID`](@ref) and [`placePI`](@ref).
"""
function loopshapingPI(P0, ωp; ϕl=0, rl=0, phasemargin=0, form::Symbol=:standard, doplot=false, Tf = nothing, F=nothing)
    issiso(P0) || throw(ArgumentError("P must be SISO"))
    if F === nothing && Tf !== nothing
        F = tf(1, [Tf^2, 2*Tf/sqrt(2), 1])
    end
    if F !== nothing
        P = P0*F
    else
        P = P0
    end
    Pw = freqresp(P, ωp)[]
    ϕp = angle(Pw)
    rp = abs.(Pw)

    if phasemargin > 0
        ϕl == 0 || @warn "Both phasemargin and ϕl provided, the provided value for ϕl will be ignored."
        ϕl = deg2rad(-180+phasemargin)
    else
        ϕl = ϕl == 0 ? ϕp : ϕl
    end
    rl = rl == 0 ? rp : rl

    kp = rl/rp*cos(ϕp-ϕl)
    ki = rl*ωp/rp*sin(ϕp-ϕl)
    C = pid(kp, ki, 0, form=:parallel)
    CF = F === nothing ? C : C*F

    fig = if doplot
        w = exp10.(LinRange(log10(ωp)-2, log10(ωp)+2, 500))
        f1 = gangoffourplot(P0,CF, w)
        f2 = nyquistplot([P0 * CF, P0], w, ylims=(-4,2), xlims=(-4,1.2), unit_circle=true, show=false, lab=["PC" "P"])
        RecipesBase.plot!([rl*cos(ϕl)], [rl*sin(ϕl)], lab="Specification point", seriestype=:scatter)
        RecipesBase.plot(f1, f2)
    else
        nothing
    end
    kp, ki = convert_pidparams_from_to(kp, ki, 0, :parallel, form)
    (; C, kp, ki, fig, CF)
end


"""
    C, kp, ki = placePI(P, ω₀, ζ; form=:standard)

Selects the parameters of a PI-controller such that the poles of 
closed loop between `P` and `C` are placed to match the poles of 
`s^2 + 2ζω₀s + ω₀^2`.

The parameters can be returned as one of several common representations 
chose by `form`, the options are
* `:standard` - ``K_p(1 + 1/(T_i s))``
* `:series` - ``K_c(1 + 1/(τ_i s))`` (equivalent to above for PI controllers)
* `:parallel` - ``K_p + K_i/s``

`C` is the returned transfer function of the controller and `params` 
is a named tuple containing the parameters.
The parameters can be accessed as `params.Kp` or `params["Kp"]` from the named tuple,
or they can be unpacked using `Kp, Ti, Td = values(params)`.

See also [`loopshapingPI`](@ref)
"""
function placePI(P::TransferFunction{<:Continuous, <:SisoRational{T}}, ω₀, ζ; form::Symbol=:standard) where T
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
    pid(kp, ki; form=:series), convert_pidparams_from_standard(kp, ki, 0, form)[1:2]...
end

placePI(sys::LTISystem, args...; kwargs...) = placePI(tf(sys), args...; kwargs...)

"""
    C, kp, ki, kd, fig, CF = loopshapingPID(P, ω; Mt = 1.3, ϕt=75, form=:standard, doplot=false, lb=-10, ub=10, Tf = 1/1000ω, F = nothing)

Selects the parameters of a PID-controller such that the Nyquist curve of the loop-transfer function ``L = PC`` at the frequency `ω` is tangent to the circle where the magnitude of ``T = PC / (1+PC)`` equals `Mt`. `ϕt` denotes the positive angle in degrees between the real axis and the tangent point.

The default values for `Mt` and `ϕt` are chosen to give a good design for processes with inertia, and may need tuning for simpler processes.

The gain of the resulting controller is generally increasing with increasing `ω` and `Mt`.

# Arguments:
- `P`: A SISO plant.
- `ω`: The specification frequency.
- `Mt`: The magnitude of the complementary sensitivity function at the specification frequency, ``|T(iω)|``.
- `ϕt`: The 
- `doplot`: If true, gang of four and Nyquist plots will be returned in `fig`.
- `lb`: log10 of lower bound for `kd`.
- `ub`: log10 of upper bound for `kd`.
- `Tf`: Time constant for second-order measurement noise filter on the form `tf(1, [Tf^2, 2*Tf/sqrt(2), 1])` to make the controller strictly proper. A practical controller typically sets this time constant slower than the default, e.g., `Tf = 1/100ω` or `Tf = 1/10ω`
- `F`: A pre-designed filter to use instead of the default second-order filter.

The parameters can be returned as one of several common representations 
chosen by `form`, the options are
* `:standard` - ``K_p(1 + 1/(T_i s) + T_ds)``
* `:series` - ``K_c(1 + 1/(τ_i s))(τ_d s + 1)``
* `:parallel` - ``K_p + K_i/s + K_d s``

See also [`loopshapingPI`](@ref), [`pidplots`](@ref), [`stabregionPID`](@ref) and [`placePI`](@ref).

# Example:
```julia
P  = tf(1, [1,0,0]) # A double integrator
Mt = 1.3  # Maximum magnitude of complementary sensitivity
ω  = 1    # Frequency at which the specification holds
C, kp, ki, kd, fig, CF = loopshapingPID(P, ω; Mt, ϕt = 75, doplot=true)
```
"""
function loopshapingPID(P0, ω; Mt = 1.3, ϕt=75, form::Symbol = :standard, doplot=false, lb=-10, ub=10, Tf = 1/1000ω, verbose=true, F=nothing)

    if F === nothing
        F = tf(1, [Tf^2, 2*Tf/sqrt(2), 1])
    end
    P = P0*F
    
    ct = -Mt^2/(Mt^2-1) # Mt center
    rt = Mt/(Mt^2-1)    # Mt radius

    specpoint = ct + rt * cis(-deg2rad(ϕt))
    rl = abs(specpoint)
    phasemargin = 180 + rad2deg(angle(specpoint))

    Pω = freqresp(P, ω)[]
    ϕp = angle(Pω)
    rp = abs.(Pω)
    dp_dω = ForwardDiff.derivative(w->freqresp(P, w)[], ω)
    ϕl = deg2rad(-180+phasemargin)

    g = rl/rp
    kp = g*cos(ϕp-ϕl)
    verbose && kp < 0 && @warn "Calculated kp is negative, try adjusting ω"
    function evalkd(_kd)
        kikd = sin(ϕp-ϕl)
        _ki = ω*(kikd + ω*_kd)
        _ki *= g
        _kd *= g

        dc_dω = complex(0, _kd + _ki/ω^2)
        Cω = kp + im*(ω*_kd - _ki/ω) # Freqresp of C
        dl_dω = Pω*dc_dω + Cω*dp_dω
        orthogonality_condition = rad2deg(angle(dl_dω)) - (90 - ϕt)
        orthogonality_condition, _ki, _kd
    end
    # Try a range of kd values to initialize bisection
    kds = exp10.(LinRange(lb, ub, 1500))
    # RecipesBase.plot(first.(evalkd.(kds))) |> display
    orths = abs.(first.(evalkd.(kds)))
    _, ind = findmin(orths)
    lb = log10(kds[max(ind-1, 1)]) # Bisect between the neighbors of the best found value in range
    ub = log10(kds[min(ind+1, end)])
    # Bisect over kd to find the root orthogonality_condition = 0
    local orthogonality_condition, ki, kd
    for i = 1:30
        midpoint = (lb+ub)/2
        kd = exp10(midpoint)
        orthogonality_condition, ki, kd = evalkd(kd)
        # @show kp, ki, kd
        orthogonality_condition
        if orthogonality_condition > 0
            lb = midpoint
        else
            ub = midpoint
        end
    end
    verbose && abs(orthogonality_condition) > 1e-5 && @warn "Bisection failed, inspect the Nyquist plot generated with doplot = true and try adjusting Mt or ϕt."
    verbose && ki < 0 && @warn "Calculated ki is negative, inspect the Nyquist plot generated with doplot = true and try adjusting ω or the angle ϕt"
    C = pid(kp, ki, kd, form=:parallel)
    any(real(p) > 0 for p in poles(C)) && @error "Calculated controller is unstable."
    kp, ki, kd = ControlSystems.convert_pidparams_from_to(kp, ki, kd, :parallel, form)
    CF = C*F
    fig = if doplot
        w = exp10.(LinRange(log10(ω)-2, log10(ω)+2, 500))
        f1 = gangoffourplot(P0,CF, w, Mt_lines=[Mt])
        f2 = nyquistplot([P0 * CF, P0], w, ylims=(-4,2), xlims=(-4,1.2), unit_circle=true, Mt_circles=[Mt], show=false, lab=["PC" "P"])
        RecipesBase.plot!([ct, real(specpoint)], [0, imag(specpoint)], lab="ϕt = $(ϕt)°", l=:dash)

        α = LinRange(0, -deg2rad(ϕt), 30)
        RecipesBase.plot!(ct .+ 0.1 .* cos.(α), 0.1 .* sin.(α), primary=false)
        RecipesBase.plot!([ct], [0], lab="T center", seriestype=:scatter, primary=false)
        RecipesBase.plot!([rl*cosd(-180+phasemargin)], [rl*sind(-180+phasemargin)], lab="Specification point", seriestype=:scatter)
        RecipesBase.plot(f1, f2)
    else
        nothing
    end
    (; C, kp, ki, kd, fig, CF)
end

"""
    Kp, Ti, Td = convert_pidparams_to_standard(param_p, param_i, param_d, form)

Convert parameters from form `form` to `:standard` form. 

The `form` can be chosen as one of the following
* `:standard` - ``K_p(1 + 1/(T_i s) + T_ds)``
* `:series` - ``K_c(1 + 1/(τ_i s))(τ_d s + 1)``
* `:parallel` - ``K_p + K_i/s + K_d s``
"""
function convert_pidparams_to_standard(param_p, param_i, param_d, form::Symbol)
    if form === :standard
        return param_p, param_i, param_d
    elseif form === :series
        return @. (
            param_p * (param_i + param_d) / param_i,
            param_i + param_d,
            param_i * param_d / (param_i + param_d)
        )
    elseif form === :parallel
        return @. (param_p, param_p / param_i, param_d / param_p)
    else
        throw(ArgumentError("form $(form) not supported."))
    end
end

"""
    param_p, param_i, param_d = convert_pidparams_from_standard(Kp, Ti, Td, form)

Convert parameters to form `form` from `:standard` form. 

The `form` can be chosen as one of the following
* `:standard` - ``K_p(1 + 1/(T_i s) + T_ds)``
* `:series` - ``K_c(1 + 1/(τ_i s))(τ_d s + 1)``
* `:parallel` - ``K_p + K_i/s + K_d s``
"""
function convert_pidparams_from_standard(Kp, Ti, Td, form::Symbol)
    if form === :standard
        return Kp, Ti, Td
    elseif form === :series
        return @. (
            (Ti - sqrt(Ti * (Ti - 4 * Td))) / 2 * Kp / Ti,
            (Ti - sqrt(Ti * (Ti - 4 * Td))) / 2,
            (Ti + sqrt(Ti * (Ti - 4 * Td))) / 2,
        )
    elseif form === :parallel
        return @. (Kp, Kp/Ti, Td*Kp)
    else
        throw(ArgumentError("form $(form) not supported."))
    end
end

"""
    convert_pidparams_from_to(kp, ki, kd, from::Symbol, to::Symbol)
"""
function convert_pidparams_from_to(kp, ki, kd, from::Symbol, to::Symbol)
    kp, ki, kd = convert_pidparams_to_standard(kp, ki, kd, from)
    convert_pidparams_from_standard(kp, ki, kd, to)
end
