export pid, pid_tf, pid_ss, pid_2dof, pid_ss_2dof, pidplots, leadlink, laglink, leadlinkat, leadlinkcurve, stabregionPID, loopshapingPI, placePI, loopshapingPID

"""
    C = pid(param_p, param_i, [param_d]; form=:standard, state_space=false, [Tf], [Ts], filter_order=2, d=1/√(2))

Calculates and returns a PID controller. 

The `form` can be chosen as one of the following (determines how the arguments `param_p, param_i, param_d` are interpreted)
* `:standard` - ``K_p(1 + 1/(T_i s) + T_d s)``
* `:series` - ``K_c(1 + 1/(τ_i s))(τ_d s + 1)``
* `:parallel` - ``K_p + K_i/s + K_d s``

If `state_space` is set to `true`, either `Kd` has to be zero
or a positive `Tf` has to be provided for creating a filter on 
the input to allow for a state-space realization. 
A balanced state-space realization is returned, unless `balance = false`.

## Filter
If `Tf` is supplied, a filter is added, the filter used is either
- `filter_order = 2` (default): ``1 / ((sT_f)^2/(4d^2) + sT_f + 1)`` in series with the controller. Note: this parametrization of the filter differs in behavior from the common parameterizaiton ``1/(s^2 + 2dws + w^2)`` as the parameters vary, the former maintains an almost fixed _bandwidth_ while `d` varies, while the latter maintains a fixed distance of the poles from the origin.
- `filter_order = 1`: ``1 / (1 + sT_f)`` applied to the derivative term only

``T_f`` can typically be chosen as ``T_i/N`` for a PI controller and ``T_d/N`` for a PID controller,
and `N` is commonly in the range 2 to 20. With a second-order filter, `d` controls the damping. `d = 1/√(2)` gives a Butterworth configuration of the poles, and `d=1` gives a critically damped filter (no overshoot). `d` above one may be used, although `d > 1` yields an increasingly over-damped filter (this parametrization does not send one pole to the origin ``d → ∞`` like the ``(ω,ζ)`` parametrization does).


## Discrete-time

For a discrete controller a positive `Ts` can be supplied.
In this case, the continuous-time controller is discretized using the Tustin method.

## Examples
```
C1 = pid(3.3, 1, 2)                             # Kd≠0 works without filter in tf form
C2 = pid(3.3, 1, 2; Tf=0.3, state_space=true)   # In statespace a filter is needed
C3 = pid(2.,  3, 0; Ts=0.4, state_space=true)   # Discrete
```

The functions `pid_tf` and `pid_ss` are also exported. They take the same parameters
and is what is actually called in `pid` based on the `state_space` parameter. See also [`pid_2dof`](@ref) for a 2DOF controller with inputs `[r; y]` and outputs `u`.
"""
function pid(param_p, param_i, param_d=zero(typeof(param_p)); form=:standard, Ts=nothing, Tf=nothing, state_space=false, balance=true, kwargs...)
    C = if state_space # Type instability? Can it be fixed easily, does it matter?
        pid_ss(param_p, param_i, param_d; form, Tf, balance, kwargs...)
    else
        pid_tf(param_p, param_i, param_d; form, Tf, kwargs...)
    end
    if Ts === nothing
        return C
    else
        param_d != 0 && Tf === nothing && throw(ArgumentError("Discretizing a continuous time PID controller without a filter is not supported. Supply a filter time constant `Tf`"))
        c2d(C, Ts, :tustin)
    end
end

@deprecate pid(; kp=0, ki=0, kd=0, series = false) pid(kp, ki, kd; form=series ? :series : :parallel)

function pid_tf(param_p, param_i, param_d=zero(typeof(param_p)); form=:standard, Tf=nothing, filter_order=2, d=1/√(2))
    Kp, Ki, Kd = convert_pidparams_to_parallel(param_p, param_i, param_d, form)
    filter_order ∈ (1,2) || throw(ArgumentError("Filter order must be 1 or 2"))
    if isnothing(Tf) || (Kd == 0 && filter_order == 1)
        if Ki == 0
            return tf([Kd, Kp], [1])
        else
            return tf([Kd, Kp, Ki], [1, 0])
        end
    else
        if Ki == 0
            if filter_order == 1
                tf([Kd*Tf + Kd, Kd], [Tf, 1])
            else
                return tf([Kd, Kp], [Tf^2/(4d^2), Tf, 1])
            end
        else
            if filter_order == 1
                return tf([Kd + Kp*Tf, Ki*Tf + Kp, Ki], [Tf, 1, 0])
            else
                return tf([Kd, Kp, Ki], [Tf^2/(4d^2), Tf, 1, 0])
            end
        end
    end
end

function pid_ss(param_p, param_i, param_d=zero(typeof(param_p)); form=:standard, Tf=nothing, balance=true, filter_order=2, d=nothing)
    Kp, Ki, Kd = convert_pidparams_to_parallel(param_p, param_i, param_d, form)
    if !isnothing(Tf)
        d42 = d === nothing ? 2.0 : 4d^2 # To avoid d = 1/sqrt(2) not yielding exactly 2
        if Ki == 0
            if filter_order == 1
                A = [-1 / Tf;;]
                B = [-Kd/Tf^2]
                C = [1.0;;]
                D = [Kd/Tf + Kp;;]
            else # 2
                A = [0 1; -d42/Tf^2 -d42/Tf]
                B = [0; 1]
                C = d42 / Tf^2 * [Kp Kd]
                D = [0.0;;]
            end
        else
            if filter_order == 1
                A = [0 0; 0 -1/Tf]
                B = [Ki; -Kd/Tf^2]
                C = [1.0 1]
                D = [Kd/Tf + Kp;;]
            else # 2
                A = [0 1 0; 0 0 1; 0 -d42/Tf^2 -d42/Tf]
                B = [0; 0; 1]
                C = d42 / Tf^2 * [Ki Kp Kd]
                D = [0.0;;]
            end
        end
    elseif Kd == 0
        if Ki != 0
            A = [0.0;;]
            B = [1.0;;]
            C = [Ki;;] # Ti == 0 would result in division by zero, but typically indicates that the user wants no integral action
            D = [Kp;;]
        else
            return ss([Kp])
        end
    else
        throw(DomainError("cannot create controller as a state space if Td != 0 without a filter. Either create the controller as a transfer function, pid(params..., state_space=false), or supply keyword argument Tf to add a filter."))
    end
    K = ss(A, B, C, D)
    balance ? first(balance_statespace(K)) : K
end

"""
    C = pid_2dof(param_p, param_i, [param_d]; form=:standard, state_space=true, N = 10, [Ts], b=1, c=0, disc=:tustin)

Calculates and returns a PID controller on 2DOF form with inputs `[r; y]` and outputs `u` where `r` is the reference signal, `y` is the measured output and `u` is the control signal.

Belowm we show two different depections of the contorller, one as a 2-input system (left) and one where the tw internal SISO systems of the controller are shown (right).
```
                                ┌──────┐                      
                             r  │      │                      
                            ───►│  Cr  ├────┐                 
r  ┌─────┐     ┌─────┐          │      │    │    ┌─────┐      
──►│     │  u  │     │ y        └──────┘    │    │     │ y    
   │  C  ├────►│  P  ├─┬─►                  +───►│  P  ├─┬───►
 ┌►│     │     │     │ │        ┌──────┐    │    │     │ │    
 │ └─────┘     └─────┘ │      y │      │    │    └─────┘ │    
 │                     │     ┌─►│  Cy  ├────┘            │    
 └─────────────────────┘     │  │      │                 │    
                             │  └──────┘                 │    
                             │                           │    
                             └───────────────────────────┘     
```

The `form` can be chosen as one of the following (determines how the arguments `param_p, param_i, param_d` are interpreted)
* `:standard` - ``K_p*(br-y + (r-y)/(T_i s) + T_d s (cr-y)/(T_f s + 1))``
* `:parallel` - ``K_p*(br-y) + K_i (r-y)/s + K_d s (cr-y)/(Tf s + 1)``

- `b` is a set-point weighting for the proportional term
- `c` is a set-point weighting for the derivative term, this defaults to 0.
- If both `b` and `c` are set to zero, the feedforward path of the controller will be strictly proper.
- `Tf` is a time constant for a filter on the derivative term, this defaults to `Td/N` where `N` is set to 10. Instead of passing `Tf` one can also pass `N` directly. The proportional term is not affected by this filter. **Please note**: this derivative filter is not the same as the one used in the `pid` function, where the filter is of second order and applied in series with the contorller, i.e., it affects all three PID terms.
- A PD controller is constructed by setting `param_i` to zero. 
- A balanced state-space realization is returned, unless `balance = false`
- If `Ts` is supplied, the controller is discretized using the method `disc` (defaults to `:tustin`).

This controller has negative feedback built in, and the closed-loop system from `r` to `y` is thus formed as
```
Cr, Cy = C[1, 1], C[1, 2]
feedback(P, Cy, pos_feedback=true)*Cr                    # Alternative 1
feedback(P, -Cy)*Cr                                      # Alternative 2
feedback(P, C, U2=2, W2=1, W1=[], pos_feedback=true) # Alternative 3, less pretty but more efficient, returns smaller realization
```
"""
function pid_2dof(args...; state_space = true, Ts = nothing, disc = :tustin, kwargs...)
    C = pid_ss_2dof(args...; kwargs...)
    Ccd = Ts === nothing ? C : c2d(C, Ts, disc)
    state_space ? Ccd : tf(Ccd)
end

function pid_ss_2dof(param_p, param_i, param_d=zero(typeof(param_p)); form=:standard, b = 1, c = 0, Tf=nothing, N=nothing, balance=true)
    kp, ki, kd = convert_pidparams_to_parallel(param_p, param_i, param_d, form)
    if kd == 0
        if ki == 0
            return ss([kp -kp])
        else
            A = [0.0;;]
            B = [ki -ki]
            C = [1.0;;] # Ti == 0 would result in division by zero, but typically indicates that the user wants no integral action
            D = [b*kp -kp]
            return ss(A, B, C, D)
        end
    end
    # On standard form we use N
    Tf !== nothing && N !== nothing && throw(ArgumentError("Cannot supply both Tf and N"))
    if Tf === nothing && N === nothing
        N = 10 # Default value
    end
    Tf = @something(Tf, kd / N)
    Tf <= 0 && throw(ArgumentError("Tf must be strictly positive"))
    if ki == 0
        A = [-(1 / Tf);;]
        B = [-kd*c/(Tf^2) kd/(Tf^2)]
        C = [1.0]
    else
        A = [0 0; 0 -(1 / Tf)]
        B = [ki -ki; -kd*c/Tf^2 kd/Tf^2]
        C = [1.0 1]
    end
    D = [kd*c/Tf+kp*b -(kd/Tf + kp)]
    K = ss(A, B, C, D)
    balance ? first(balance_statespace(K)) : K
end


"""
    pidplots(P, args...; params_p, params_i, params_d=0, form=:standard, ω=0, grid=false, kwargs...)

Display the relevant plots related to closing the loop around process `P` with a PID controller supplied in `params`
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
        bodeplot(Cs, ω; lab=repeat(labels, inner=(1,2)), title="Controller bode plot", kwargs...) |> display
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
    leadlink(b, N, K=1; [Ts])

Returns a phase advancing link, the top of the phase curve is located at `ω = b√(N)` where the link amplification is `K√(N)` The bode curve will go from `K`, bend up at `b` and level out at `KN` for frequencies > `bN`

The phase advance at `ω = b√(N)` can be plotted as a function of `N` with `leadlinkcurve()`

Values of `N < 1` will give a phase retarding link.

```math
KN \\dfrac{s + b}{s + bN} = K \\dfrac{1 + s/b}{1 + s/(bN)}
```

See also `leadlinkat` `laglink`
"""
function leadlink(b, N, K=1; Ts=nothing)
    Ts !== nothing && (Ts ≥ 0 || throw(ArgumentError("Negative `Ts` is not supported.")))
    N > 1 || @warn "N should be ≥ 1 for the link to be phase advancing."
    numerator = [1/b, 1]
    denominator = [1/(b*N), 1]
    gain = K
    G = tf(gain*numerator,denominator)
    return  isnothing(Ts) ? G : c2d(G,Ts)
end

"""
    leadlinkat(ω, N, K=1; [Ts])

Returns a phase advancing link, the top of the phase curve is located at `ω` where the link amplification is `K√(N)` The bode curve will go from `K`, bend up at `ω/√(N)` and level out at `KN` for frequencies > `ω√(N)`

The phase advance at `ω` can be plotted as a function of `N` with `leadlinkcurve()`

Values of `N < 1` will give a phase retarding link.

See also `leadlink` `laglink`
"""
function leadlinkat(ω, N, K=1; Ts=nothing)
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
    K   = convert_pidparams_from_parallel.(kp, ki, kd, form)
    kp, ki = getindex.(K, 1), getindex.(K, 2)
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
    K       = convert_pidparams_from_parallel.(kp, ki, kd, form)
    kp, ki  = getindex.(K, 1), getindex.(K, 2)
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
* `:standard` - ``K_p(1 + 1/(T_i s) + T_d s)``
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
    kp, ki = convert_pidparams_from_parallel(kp, ki, 0, form)
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
- `ϕt`: The positive angle in degrees between the real axis and the tangent point.
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
    iscontinuous(P0) || throw(ArgumentError("Discrete-time system models are not supported, convert to continuous time with d2c(P)"))
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
    kp, ki, kd = convert_pidparams_from_parallel(kp, ki, kd, form)
    CF = C*F
    fig = if doplot
        w = exp10.(LinRange(log10(ω)-2, log10(ω)+2, 1000))
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
        return (param_p, param_i, param_d)
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
    Kp, Ti, Td = convert_pidparams_to_parallel(param_p, param_i, param_d, form)

Convert parameters from form `form` to `:parallel` form.

The `form` can be chosen as one of the following
* `:standard` - ``K_p(1 + 1/(T_i s) + T_d s)``
* `:series` - ``K_c(1 + 1/(τ_i s))(τ_d s + 1)``
* `:parallel` - ``K_p + K_i/s + K_d s``
"""
function convert_pidparams_to_parallel(param_p, param_i, param_d, form::Symbol)
    if form === :parallel
        return (param_p, param_i, param_d)
    elseif form === :series
        # param_i = 0 would result in division by zero, but typically indicates that the user wants no integral action
        param_i == 0 && return (param_p, 0, param_p * param_d)
        return (param_p * (param_i + param_d) / param_i,
                param_p / param_i,
                param_p * param_d)
    elseif form === :standard
        param_i == 0 && return (param_p, 0, param_p * param_d)
        return (param_p, param_p / param_i, param_p * param_d)
    else
        throw(ArgumentError("form $(form) not supported."))
    end
end

"""
    param_p, param_i, param_d = convert_pidparams_from_standard(Kp, Ti, Td, form)

Convert parameters to form `form` from `:standard` form. 

The `form` can be chosen as one of the following
* `:standard` - ``K_p(1 + 1/(T_i s) + T_d s)``
* `:series` - ``K_c(1 + 1/(τ_i s))(τ_d s + 1)``
* `:parallel` - ``K_p + K_i/s + K_d s``
"""
function convert_pidparams_from_standard(Kp, Ti, Td, form::Symbol)
    if form === :standard
        return (Kp, Ti, Td)
    elseif form === :series
        Δ = Ti * (Ti - 4 * Td)
        Δ < 0 && throw(DomainError("The condition Ti^2 ≥ 4Td*Ti is not satisfied: the PID parameters cannot be converted to series form"))
        sqrtΔ = sqrt(Δ)
        return ((Ti - sqrtΔ) / 2 * Kp / Ti,
                (Ti - sqrtΔ) / 2,
                (Ti + sqrtΔ) / 2)
    elseif form === :parallel
        return (Kp, Kp/Ti, Td*Kp)
    else
        throw(ArgumentError("form $(form) not supported."))
    end
end


"""
    Kp, Ti, Td = convert_pidparams_from_parallel(Kp, Ki, Kd, to_form)

Convert parameters from form `:parallel` to form `to_form`.

The `form` can be chosen as one of the following
* `:standard` - ``K_p(1 + 1/(T_i s) + T_d s)``
* `:series` - ``K_c(1 + 1/(τ_i s))(τ_d s + 1)``
* `:parallel` - ``K_p + K_i/s + K_d s``
"""
function convert_pidparams_from_parallel(Kp, Ki, Kd, to::Symbol)
    if to === :parallel
        return (Kp, Ki, Kd)
    elseif to === :series
        Ki == 0 && return (Kp, 0, Kp*Kd)
        Δ = Kp^2-4Ki*Kd
        Δ < 0 &&
            throw(DomainError("The condition Kp^2 ≥ 4Ki*Kd is not satisfied: the PID parameters cannot be converted to series form"))
        sqrtΔ = sqrt(Δ)
        return ((Kp - sqrtΔ)/2, (Kp - sqrtΔ)/(2Ki), (Kp + sqrtΔ)/(2Ki))
    elseif to === :standard
        Kp == 0 && throw(DomainError("Cannot convert to standard form when Kp=0"))
        Ki == 0 && return (Kp, Inf, Kd / Kp)
        return (Kp, Kp / Ki, Kd / Kp)
    else
        throw(ArgumentError("form $(form) not supported."))
    end
end


"""
    convert_pidparams_from_to(kp, ki, kd, from::Symbol, to::Symbol)
"""
function convert_pidparams_from_to(kp, ki, kd, from::Symbol, to::Symbol)
    Kp, Ki, Kd = convert_pidparams_to_parallel(kp, ki, kd, from)
    convert_pidparams_from_parallel(Kp, Ki, Kd, to)
end
