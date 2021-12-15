"""
    poles(sys)

Compute the poles of system `sys`."""
poles(sys::AbstractStateSpace) = eigvalsnosort(sys.A)
poles(sys::SisoTf) = error("pole is not implemented for type $(typeof(sys))")

# Seems to have a lot of rounding problems if we run the full thing with sisorational,
# converting to zpk before works better in the cases I have tested.
poles(sys::TransferFunction) = poles(zpk(sys))

function poles(sys::TransferFunction{<:TimeEvolution,SisoZpk{T,TR}}) where {T, TR}
    # With right TR, this code works for any SisoTf

    # Calculate least common denominator of the minors,
    # i.e. something like least common multiple of the pole-polynomials
    individualpoles = [map(poles, sys.matrix)...;]
    lcmpoles = TR[]
    for poles = minorpoles(sys.matrix)
        # Poles have to be equal to existing poles for the individual transfer functions and this
        # calculation probably is more precise than the full. Seems to work better at least.
        for i = 1:length(poles)
            idx = argmin(map(abs, individualpoles .- poles[i]))
            poles[i] = individualpoles[idx]
        end
        for pole = lcmpoles
            idx = findfirst(poles .≈ pole)
            if idx != nothing
                deleteat!(poles, idx)
            end
        end
        append!(lcmpoles, poles)
    end

    return lcmpoles
end

"""
    minorpoles(sys)

Compute the poles of all minors of the system."""
# TODO: Improve implementation, should be more efficient ways.
# Calculates the same minors several times in some cases.
function minorpoles(sys::Matrix{SisoZpk{T, TR}}) where {T, TR}
    minors = Array{TR,1}[]
    ny, nu = size(sys)
    if ny == nu == 1
        push!(minors, poles(sys[1, 1]))
    elseif ny == nu
        push!(minors, poles(det(sys)))
        for i = 1:ny
            for j = 1:nu
                newmat = sys[1:end .!=i, 1:end .!= j]
                append!(minors, minorpoles(newmat))
            end
        end
    elseif ny < nu
        for i = 1:nu
            newmat = sys[1:end, 1:end .!= i]
            append!(minors, minorpoles(newmat))
        end
    else
        for i = 1:ny
            newmat = sys[1:end .!= i, 1:end]
            append!(minors, minorpoles(newmat))
        end
    end
    return minors
end

"""
    det(sys)

Compute the determinant of the Matrix `sys` of SisoTf systems, returns a SisoTf system."""
# TODO: improve this implementation, should be more efficient ones
function det(sys::Matrix{S}) where {S<:SisoZpk}
    ny, nu = size(sys)
    ny == nu || throw(ArgumentError("sys matrix is not square"))
    if ny == 1
        return sys[1, 1]
    end
    tot = zero(S)
    sign = -1
    for i = 1:ny
        sign = -sign
        tot += sign * sys[i, 1] * det(sys[1:end .!= i, 2:end])
    end
    return tot
end

"""
    dcgain(sys, ϵ=0)

Compute the dcgain of system `sys`.

equal to G(0) for continuous-time systems and G(1) for discrete-time systems.

`ϵ` can be provided to evaluate the dcgain with a small perturbation into
the stability region of the complex plane.
"""
function dcgain(sys::LTISystem, ϵ=0)
    return iscontinuous(sys) ? evalfr(sys, -ϵ) : evalfr(sys, exp(-ϵ*sys.Ts))
end

"""
    markovparam(sys, n)

Compute the `n`th markov parameter of state-space system `sys`. This is defined
as the following:

`h(0) = D`

`h(n) = C*A^(n-1)*B`"""
function markovparam(sys::AbstractStateSpace, n::Integer)
    n < 0 && error("n must be >= 0")
    return n == 0 ? sys.D : sys.C * sys.A^(n-1) * sys.B
end

"""
    z, p, k = zpkdata(sys)

Compute the zeros, poles, and gains of system `sys`.

### Returns
`z` : Matrix{Vector{ComplexF64}}, (ny x nu)

`p` : Matrix{Vector{ComplexF64}}, (ny x nu)

`k` : Matrix{Float64}, (ny x nu)"""
function zpkdata(sys::LTISystem)
    G = convert(TransferFunction{typeof(timeevol(sys)),SisoZpk}, sys)

    zs = map(x -> x.z, G.matrix)
    ps = map(x -> x.p, G.matrix)
    ks = map(x -> x.k, G.matrix)

    return zs, ps, ks
end

"""
    Wn, zeta, ps = damp(sys)

Compute the natural frequencies, `Wn`, and damping ratios, `zeta`, of the
poles, `ps`, of `sys`"""
function damp(sys::LTISystem)
    ps = poles(sys)
    if isdiscrete(sys)
        ps = log.(complex.(ps))/sys.Ts
    end
    Wn = abs.(ps)
    order = sortperm(Wn; by=z->(abs(z), real(z), imag(z)))
    Wn = Wn[order]
    ps = ps[order]
    ζ = -cos.(angle.(ps))
    return Wn, ζ, ps
end

"""
    dampreport(sys)

Display a report of the poles, damping ratio, natural frequency, and time
constant of the system `sys`"""
function dampreport(io::IO, sys::LTISystem)
    Wn, zeta, ps = damp(sys)
    t_const = 1 ./ (Wn.*zeta)
    header =
    ("|        Pole        |   Damping     |   Frequency   |   Frequency   | Time Constant |\n"*
     "|                    |    Ratio      |   (rad/sec)   |     (Hz)      |     (sec)     |\n"*
     "+--------------------+---------------+---------------+---------------+---------------+")
    println(io, header)
    if all(isreal, ps)
        for i=eachindex(ps)
            p, z, w, t = ps[i], zeta[i], Wn[i], t_const[i]
            Printf.@printf(io, "| %-+18.3g |  %-13.3g|  %-13.3g|  %-13.3g|  %-13.3g|\n", real(p), z, w, w/(2π), t)
        end
    elseif numeric_type(sys) <: Real # real-coeff system with complex conj. poles
        for i=eachindex(ps)
            p, z, w, t = ps[i], zeta[i], Wn[i], t_const[i]
            imag(p) < 0 && (continue) # use only the positive complex pole to print with the ± operator
            if imag(p) == 0 # no ± operator for real pole
                Printf.@printf(io, "| %-+18.3g |  %-13.3g|  %-13.3g|  %-13.3g|  %-13.3g|\n", real(p), z, w, w/(2π), t)
            else
                Printf.@printf(io, "| %-+7.3g ± %6.3gim |  %-13.3g|  %-13.3g|  %-13.3g|  %-13.3g|\n", real(p), imag(p), z, w, w/(2π), t)
            end
        end
    else # complex-coeff system
        for i=eachindex(ps)
            p, z, w, t = ps[i], zeta[i], Wn[i], t_const[i]
            Printf.@printf(io, "| %-+7.3g  %+7.3gim |  %-13.3g|  %-13.3g|  %-13.3g|  %-13.3g|\n", real(p), imag(p), z, w, w/(2π), t)
        end
    end
end
dampreport(sys::LTISystem) = dampreport(stdout, sys)


"""
    tzeros(sys)

Compute the invariant zeros of the system `sys`. If `sys` is a minimal
realization, these are also the transmission zeros."""
function tzeros(sys::TransferFunction)
    if issiso(sys)
        return tzeros(sys.matrix[1,1])
    else
        return tzeros(ss(sys))
    end
end

# Implements the algorithm described in:
# Emami-Naeini, A. and P. Van Dooren, "Computation of Zeros of Linear
# Multivariable Systems," Automatica, 18 (1982), pp. 415–430.
#
# Note that this returns either Vector{ComplexF32} or Vector{Float64}
tzeros(sys::AbstractStateSpace) = tzeros(sys.A, sys.B, sys.C, sys.D)
# Make sure everything is BlasFloat
function tzeros(A::AbstractMatrix, B::AbstractMatrix, C::AbstractMatrix, D::AbstractMatrix)
    T = promote_type(eltype(A), eltype(B), eltype(C), eltype(D))
    A2, B2, C2, D2, _ = promote(A,B,C,D, fill(zero(T)/one(T),0,0)) # If Int, we get Float64
    tzeros(A2, B2, C2, D2)
end
function tzeros(A::AbstractMatrix{T}, B::AbstractMatrix{T}, C::AbstractMatrix{T}, D::AbstractMatrix{T}) where {T <: Union{AbstractFloat,Complex{<:AbstractFloat}}#= For eps(T) =#}
    # Balance the system
    A, B, C = balance_statespace(A, B, C)

    # Compute a good tolerance
    meps = 10*eps(real(T))*norm([A B; C D])

    # Step 1:
    A_r, B_r, C_r, D_r = reduce_sys(A, B, C, D, meps)

    # Step 2: (conjugate transpose should be avoided since single complex zeros get conjugated)
    A_rc, B_rc, C_rc, D_rc = reduce_sys(transpose(A_r), transpose(C_r), transpose(B_r), transpose(D_r), meps)
    if isempty(A)   return complex(T)[]    end

    # Step 3:
    # Compress cols of [C D] to [0 Df]
    mat = [C_rc D_rc]
    # To ensure type-stability, we have to annote the type here, as qrfact
    # returns many different types.
    W = qr(mat').Q
    W = reverse(W, dims=2)
    mat = mat*W
    if fastrank(mat', meps) > 0
        nf = size(A_rc, 1)
        m = size(D_rc, 2)
        Af = ([A_rc B_rc] * W)[1:nf, 1:nf]
        Bf = ([Matrix{T}(I, nf, nf) zeros(nf, m)] * W)[1:nf, 1:nf]
        zs = eigvalsnosort(Af, Bf)
        _fix_conjugate_pairs!(zs) # Generalized eigvals does not return exact conj. pairs
    else
        zs = complex(T)[]
    end
    return zs
end


"""
Implements REDUCE in the Emami-Naeini & Van Dooren paper. Returns transformed
A, B, C, D matrices. These are empty if there are no zeros.
"""
function reduce_sys(A::AbstractMatrix, B::AbstractMatrix, C::AbstractMatrix, D::AbstractMatrix, meps::AbstractFloat)
    T = promote_type(eltype(A), eltype(B), eltype(C), eltype(D))
    Cbar, Dbar = C, D
    if isempty(A)
        return A, B, C, D
    end
    while true
        # Compress rows of D
        U = qr(D).Q
        D = U'D
        C = U'C
        sigma = fastrank(D, meps)
        Cbar = C[1:sigma, :]
        Dbar = D[1:sigma, :]
        Ctilde = C[(1 + sigma):end, :]
        if sigma == size(D, 1)
            break
        end

        # Compress columns of Ctilde
        V = qr(Ctilde').Q
        V = reverse(V,dims=2)
        Sj = Ctilde*V
        rho = fastrank(Sj', meps)
        nu = size(Sj, 2) - rho

        if rho == 0
            break
        elseif nu == 0
            # System has no zeros, return empty matrices
            A = B = Cbar = Dbar = Matrix{T}(undef, 0,0)
            break
        end
        # Update System
        n, m = size(B)
        Vm = [V zeros(T, n, m); zeros(T, m, n) Matrix{T}(I, m, m)]
        if sigma > 0
            M = [A B; Cbar Dbar]
            Vs = [V' zeros(T, n, sigma) ; zeros(T, sigma, n) Matrix{T}(I, sigma, sigma)]
        else
            M = [A B]
            Vs = V'
        end
        sigma, rho, nu
        M = Vs * M * Vm
        A = M[1:nu, 1:nu]
        B = M[1:nu, (nu + rho + 1):end]
        C = M[(nu + 1):end, 1:nu]
        D = M[(nu + 1):end,  (nu + rho + 1):end]
    end
    return A, B, Cbar, Dbar
end

# Determine the number of non-zero rows, with meps as a tolerance. For an
# upper-triangular matrix, this is a good proxy for determining the row-rank.
function fastrank(A::AbstractMatrix, meps::Real)
    n, m = size(A)
    if n*m == 0     return 0    end
    norms = Vector{real(eltype(A))}(undef, n)
    for i = 1:n
        norms[i] = norm(A[i, :])
    end
    mrank = sum(norms .> meps)
    return mrank
end

"""
`ωgₘ, gₘ, ωϕₘ, ϕₘ = margin{S<:Real}(sys::LTISystem, w::AbstractVector{S}; full=false, allMargins=false)`

returns frequencies for gain margins, gain margins, frequencies for phase margins, phase margins

If `!allMargins`, return only the smallest margin

If `full` return also `fullPhase`

See also [`diskmargin`](@ref)
"""
function margin(sys::LTISystem, w::AbstractVector{<:Real}; full=false, allMargins=false)
    ny, nu = size(sys)

    if allMargins
        wgm         = Array{Array{numeric_type(sys),1}}(undef, ny,nu)
        gm          = Array{Array{numeric_type(sys),1}}(undef, ny,nu)
        wpm         = Array{Array{numeric_type(sys),1}}(undef, ny,nu)
        pm          = Array{Array{numeric_type(sys),1}}(undef, ny,nu)
        fullPhase   = Array{Array{numeric_type(sys),1}}(undef, ny,nu)
    else
        wgm         = Array{numeric_type(sys),2}(undef, ny, nu)
        gm          = Array{numeric_type(sys),2}(undef, ny, nu)
        wpm         = Array{numeric_type(sys),2}(undef, ny, nu)
        pm          = Array{numeric_type(sys),2}(undef, ny, nu)
        fullPhase   = Array{numeric_type(sys),2}(undef, ny, nu)
    end
    for j=1:nu
        for i=1:ny
            wgm[i,j], gm[i,j], wpm[i,j], pm[i,j], fullPhase[i,j] = sisomargin(sys[i,j], w, full=true, allMargins=allMargins)
        end
    end
    if full
        wgm, gm, wpm, pm, fullPhase
    else
        wgm, gm, wpm, pm
    end
end

"""
`ωgₘ, gₘ, ωϕₘ, ϕₘ = sisomargin{S<:Real}(sys::LTISystem, w::AbstractVector{S}; full=false, allMargins=false)`

returns frequencies for gain margins, gain margins, frequencies for phase margins, phase margins. See also [`diskmargin`](@ref).
"""
function sisomargin(sys::LTISystem, w::AbstractVector{<:Real}; full=false, allMargins=false)
    ny, nu = size(sys)
    if ny !=1 || nu != 1
        error("System must be SISO, use `margin` instead")
    end
    mag, phase, w = bode(sys, w)
    wgm, = _allPhaseCrossings(w, phase)
    gm = similar(wgm)
    for i = eachindex(wgm)
        gm[i] = 1 ./ abs(freqresp(sys,[wgm[i]])[1][1])
    end
    wpm, fi = _allGainCrossings(w, mag)
    pm = similar(wpm)
    for i = eachindex(wpm)
        pm[i] = mod(rad2deg(angle(freqresp(sys,[wpm[i]])[1][1])),360)-180
    end
    if !allMargins #Only output the smallest margins
        gm, idx = findmin([gm;Inf])
        wgm = [wgm;NaN][idx]
        fi = [fi;NaN][idx]
        pm, idx = findmin([abs.(pm);Inf])
        wpm = [wpm;NaN][idx]
        if full
            if !isnan(fi) #fi may be NaN, fullPhase is a scalar
                fullPhase = interpolate(fi, phase)
            else
                fullPhase = NaN
            end
        end
    else
        if full #We know that all values are defined and fullPhase is a vector
            fullPhase = interpolate(fi, phase)
        end
    end
    if full
        wgm, gm, wpm, pm, fullPhase
    else
        wgm, gm, wpm, pm
    end
end
margin(system::LTISystem; kwargs...) =
margin(system, _default_freq_vector(system, Val{:bode}()); kwargs...)
#margin(sys::LTISystem, args...) = margin(LTISystem[sys], args...)

# Interpolate the values in "list" given the floating point "index" fi
function interpolate(fi, list)
    fif = floor.(Integer, fi)
    fic = ceil.(Integer, fi)
    list[fif]+mod.(fi,1).*(list[fic]-list[fif])
end

function _allGainCrossings(w, mag)
    _findCrossings(w,mag.>1,mag.-1)
end

function _allPhaseCrossings(w, phase)
    #Calculate numer of times real axis is crossed on negative side
    n = fld.(phase.+180,360) #Nbr of crossed
    ph = mod.(phase,360) .- 180 #Residual
    _findCrossings(w, n, ph)
end

function _findCrossings(w, n, res)
    wcross = Vector{eltype(w)}()
    tcross = Vector{eltype(w)}()
    for i in 1:(length(w)-1)
        if res[i] == 0
            wcross = [wcross; w[i]]
            tcross = [tcross; i]
        elseif n[i] != n[i+1]
            #Interpolate to approximate crossing
            t = res[i]/(res[i]-res[i+1])
            tcross = [tcross; i+t]
            wt = w[i]+t*(w[i+1]-w[i])
            wcross = [wcross; wt]
        end
    end
    if res[end] == 0 #Special case if multiple points
        wcross = [wcross; w[end]]
        tcross = [tcross; length(w)]
    end
    wcross, tcross
end

"""
    dₘ = delaymargin(G::LTISystem)

Only supports SISO systems

See also [`diskmargin`](@ref)
"""
function delaymargin(G::LTISystem)
    # Phase margin in radians divided by cross-over frequency in rad/s.
    if G.nu + G.ny > 2
        error("delaymargin only supports SISO systems")
    end
    m     = margin(G,allMargins=true)
    ϕₘ, i = findmin(m[4])
    ϕₘ   *= π/180
    ωϕₘ   = m[3][i]
    dₘ    = ϕₘ/ωϕₘ
    if isdiscrete(G)
        dₘ /= G.Ts # Give delay margin in number of sample times, as matlab does
    end
    dₘ
end

"""
    Diskmargin

The notation follows "An Introduction to Disk Margins", Peter Seiler, Andrew Packard, and Pascal Gahinet

# Fields:
`α`: The disk margin
`ω0`: The worst-case frequency
`f0`: The destabilizing perturbation `f0` is a complex number with simultaneous gain and phase variation. This critical perturbation causes an instability with closed-loop pole on the imaginary axis at the critical frequency ω0 
`δ0`: The uncertain element generating f0.
`γmin`: The lower real-axis intercept of the disk (classical lower gain margin).
`γmax`: The upper real-axis intercept of the disk (classical upper gain margin).
`ϕm`: is the classical phase margin.
`σ`: The skew parameter that was used to calculate the margin
"""
struct Diskmargin
    α::Float64
    ω0::Float64
    f0::ComplexF64
    δ0::ComplexF64
    γmin::Float64
    γmax::Float64
    σ::Float64
    ϕm::Float64
    L
end

function Base.show(io::IO, dm::Diskmargin)
    println(io, "Disk margin with:")
    println(io, "Margin: ", dm.α)
    println(io, "Frequency: ", dm.ω0)
    println(io, "Gain margins: [$(dm.γmin), $(dm.γmax)]")
    println(io, "Phase margin: ", dm.ϕm)
    println(io, "Skew: ", dm.σ)
    println(io, "Worst-case perturbation: ", dm.f0)
end

struct Disk
    γmin::Float64
    γmax::Float64
    c::Float64
    r::Float64
    ϕm::Float64
end

center_radius(γmin, γmax) = 1/2 * (γmax + γmin), 1/2 * (γmax - γmin)

function Disk(γmin, γmax)
    c, r = center_radius(γmin, γmax)
    Disk(γmin, γmax, c, r)
end

function Disk(γmin, γmax, c, r)
    if !isfinite(γmax)
        ϕm = 90
    else
        ϕm = (1 + γmin*γmax) / (γmin + γmax)
        ϕm = ϕm >= 1 ? Inf : rad2deg(acos(ϕm))
    end
    Disk(γmin, γmax, c, r, ϕm)
end

function Disk(; α, σ)
    if α >= 2/abs(1 + σ) # the case α ≥ 2 / |1+σ| can be used to model situations where the gain can vary substantially or the phase is essentially unknown.
        # In this case, we have an "inverted circle", but we return a more conservative region corresponding to an infinite circle extending to the right
        γmin = max(0, (2 + α*(1-σ)) / (2 - α*(1+σ))) # expression for γmax in the normal case. 
        γmax = Inf
    else
        γmin = max(0, (2 - α*(1-σ)) / (2 + α*(1+σ)))
        γmax = (2 + α*(1-σ)) / (2 - α*(1+σ))
    end
    Disk(γmin, γmax)
end

Disk(dm::Diskmargin) = Disk(dm.γmin, dm.γmax)

"""
    diskmargin(L, σ = 0)
    diskmargin(L, σ::Real, ω)

Calculate the disk margin of LTI system `L`.

The implementation and notation follows
"An Introduction to Disk Margins", Peter Seiler, Andrew Packard, and Pascal Gahinet

The margins are aviable as fields of the returned objects, see [`Diskmargin`](@ref).

# Arguments:
- `L`: A loop-transfer function.
- `σ`: If little is known about the distribution of gain variations then σ = 0
is a reasonable choice as it allows for a gain increase or decrease by the same relative amount.
The choice σ < 0 is justified if the gain can decrease by a larger factor than it can increase.
Similarly, the choice σ > 0 is justified when the gain can increase by a larger factor than it can
decrease.
- `kwargs`: Are sent to the [`hinfnorm`](@ref) calculation
- `ω`: If a vector of frequencies is supplied, the frequency-dependent disk margin will be computed, see example below.

# Example: 
```
L = tf(25, [1,10,10,10])
dm = diskmargin(L, 0)
plot(dm) # Plot the disk margin to illustrate maximum allowed simultaneous gain and phase variations.

nyquistplot(L)
plot!(dm, nyquist=true) # plot a nyquist exclusion disk. The Nyquist curve will be tangent to this disk at `dm.ω0`
nyquistplot!(dm.f0*L) # If we perturb the system with the worst-case perturbation `f0`, the curve will pass through the critical point -1.

## Frequency-dependent margin
w = exp10.(LinRange(-2, 2, 500))
dms = diskmargin(L, 0, w)
plot(w, dms)
```
"""
function diskmargin(L, σ::Real=0; kwargs...)
    issiso(L) || error("MIMO not yet supported in diskmargin.")
    S̄ = 1/(1 + L) + (σ-1)/2
    n,ω0 = hinfnorm(S̄; kwargs...)
    diskmargin(L, σ, ω0)
end

diskmargin(L, σ::Real, ω::AbstractArray) = map(w->diskmargin(L, σ, w), ω)

function diskmargin(L, σ::Real, ω0::Real)
    issiso(L) || error("MIMO not yet supported in diskmargin.") # Calculation of structured singular values required, determine if det(I-MΔ) can be 0
    S̄ = 1/(1 + L) + (σ-1)/2
    freq = isdiscrete(L) ? cis(ω0*L.Ts) : complex(0, ω0)
    Sω = S̄(freq)[]
    αmax = 1/abs(Sω)
    δ0 = inv(Sω)
    dp = Disk(; α = αmax, σ)
    if δ0 == 2/(σ+1)  # S = 1, L = 0
        Diskmargin(αmax, ω0, Inf, δ0, 0, Inf, σ, dp.ϕm, L)
    else
        f0 = (2 + δ0*(1-σ)) / (2 - δ0*(1+σ))
        Diskmargin(αmax, ω0, f0, δ0, dp.γmin, dp.γmax, σ, dp.ϕm, L)
    end
end


function robust_minreal(G, args...; kwargs...)
    try
        return minreal(G, args...; kwargs...)
    catch
        return G
    end
end

"""
    S, PS, CS, T = gangoffour(P, C; minimal=true)
    gangoffour(P::AbstractVector, C::AbstractVector; minimal=true)

Given a transfer function describing the plant `P` and a transfer function describing the controller `C`, computes the four transfer functions in the Gang-of-Four.

- `S = 1/(1+PC)` Sensitivity function
- `PS = P/(1+PC)` Load disturbance to measurement signal
- `CS = C/(1+PC)` Measurement noise to control signal
- `T = PC/(1+PC)` Complementary sensitivity function

Only supports SISO systems
"""
function gangoffour(P::LTISystem, C::LTISystem; minimal=true)
    if !issiso(P) || !issiso(C)
        error("gangoffour only supports SISO systems")
    end
    minfun = minimal ? robust_minreal : identity
    S = feedback(1, P*C)    |> minfun
    PS = feedback(P, C)     |> minfun
    CS = feedback(C, P)     |> minfun
    T = feedback(P*C, 1)    |> minfun
    return S, PS, CS, T
end

"""
    S, PS, CS, T, RY, RU, RE = gangofseven(P,C,F)

Given transfer functions describing the Plant `P`, the controller `C` and a feed forward block `F`,
computes the four transfer functions in the Gang-of-Four and the transferfunctions corresponding to the feed forward.

`S = 1/(1+PC)` Sensitivity function

`PS = P/(1+PC)`

`CS = C/(1+PC)`

`T = PC/(1+PC)` Complementary sensitivity function

`RY = PCF/(1+PC)`

`RU = CF/(1+P*C)`

`RE = F/(1+P*C)`

Only supports SISO systems"""
function gangofseven(P::TransferFunction, C::TransferFunction, F::TransferFunction)
    if !issiso(P) || !issiso(C) || !issiso(F)
        error("gangofseven only supports SISO systems")
    end
    S, PS, CS, T = gangoffour(P,C)
    RY = T*F
    RU = CS*F
    RE = S*F
    return S, PS, CS, T, RY, RU, RE
end
