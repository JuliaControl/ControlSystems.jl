"""
    poles(sys)

Compute the poles of system `sys`.

Note: Poles with multiplicity `n > 1` may suffer numerical inaccuracies on the order `eps(numeric_type(sys))^(1/n)`, i.e., a double pole in a system with `Float64` coefficients may be computed with an error of about `√(eps(Float64)) ≈ 1.5e-8`.

To compute the poles of a system with non-BLAS floats, such as `BigFloat`, install and load the package `GenericSchur.jl` before calling `poles`.
"""
poles(sys::AbstractStateSpace) = eigvalsnosort(sys.A)

# Seems to have a lot of rounding problems if we run the full thing with sisorational,
# converting to zpk before works better in the cases I have tested.
function poles(sys::TransferFunction{<:TimeEvolution,SisoRational{T}}; kwargs...) where {T}
    n,d = numpoly(sys), denpoly(sys)
    sort!(filter!(isfinite, MatrixPencils.rmpoles(n, d; kwargs...)[1]), by=abs)
end

function poles(sys::TransferFunction{<:TimeEvolution,SisoZpk{T,TR}}) where {T, TR}
    # With right TR, this code works for any SisoTf

    # Calculate least common denominator of the minors,
    # i.e. something like least common multiple of the pole-polynomials
    individualpoles = [map(poles, sys.matrix)...;]
    lcmpoles = TR[]
    for poles = minorpoles(sys.matrix)
        # Poles have to be equal to existing poles for the individual transfer functions and this
        # calculation probably is more precise than the full. Seems to work better at least.
        for i = eachindex(poles)
            idx = argmin(map(abs, individualpoles .- poles[i]))
            poles[i] = individualpoles[idx]
        end
        for pole = lcmpoles
            idx = findfirst(poles .≈ pole)
            if idx !== nothing
                deleteat!(poles, idx)
            end
        end
        append!(lcmpoles, poles)
    end

    return lcmpoles
end


function count_eigval_multiplicity(p, location, e=eps(maximum(abs, p)))
    n = length(p)
    n == 0 && return 0
    for i = 1:n
        # if we count i poles within the circle assuming i integrators, we return i
        if count(p->abs(p-location) < (e^(1/i)), p) == i
            return i
        end
    end
    0
end

"""
    count_integrators(P)

Count the number of poles in the origin by finding the first value of `n` for which the number of poles within a circle of radius `eps(maximum(abs, p))^(1/n)` around the origin (1 in discrete time) equals `n`.

See also [`integrator_excess`](@ref).
"""
function count_integrators(P::LTISystem)
    p = poles(P)
    location = iscontinuous(P) ? 0 : 1
    count_eigval_multiplicity(p, location)
end

"""
    integrator_excess(P)

Count the number of integrators in the system by finding the difference between the number of poles in the origin and the number of zeros in the origin. If the number of zeros in the origin is greater than the number of poles in the origin, the count is negative.

For discrete-time systems, the origin ``s = 0`` is replaced by the point ``z = 1``.
"""
function integrator_excess(P::LTISystem)
    p = poles(P)
    z = tzeros(P)
    location = iscontinuous(P) ? 0 : 1
    count_eigval_multiplicity(p, location) - count_eigval_multiplicity(z, location)
end


# TODO: Improve implementation, should be more efficient ways.
# Calculates the same minors several times in some cases.
"""
    minorpoles(sys)

Compute the poles of all minors of the system."""
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

# TODO: improve this implementation, should be more efficient ones
"""
    det(sys)

Compute the determinant of the Matrix `sys` of SisoTf systems, returns a SisoTf system."""
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
dcgain(G::Union{UniformScaling, Number, AbstractMatrix}) = G

"""
    markovparam(sys, n)

Compute the `n`th markov parameter of discrete-time state-space system `sys`. This is defined
as the following:

`h(0) = D`

`h(n) = C*A^(n-1)*B`"""
function markovparam(sys::AbstractStateSpace{<:Discrete}, n::Integer)
    n < 0 && error("n must be >= 0")
    return n == 0 ? sys.D : sys.C * sys.A^(n-1) * sys.B
end

"""
    z, p, k = zpkdata(sys)

Compute the zeros, poles, and gains of system `sys`.

### Returns
- `z` : Matrix{Vector{ComplexF64}}, (ny × nu)
- `p` : Matrix{Vector{ComplexF64}}, (ny × nu)
- `k` : Matrix{Float64}, (ny × nu)"""
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
        ps = @. log(complex(ps))/sys.Ts
    end
    Wn = abs.(ps)
    order = sortperm(Wn; by=z->(abs(z), real(z), imag(z)))
    Wn = Wn[order]
    ps = ps[order]
    ζ = @. -cos(angle(ps))
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
    tzeros(sys::AbstractStateSpace; extra=Val(false))

Compute the invariant zeros of the system `sys`. If `sys` is a minimal
realization, these are also the transmission zeros.

If `sys` is a state-space system the function has additional keyword arguments, see [`?ControlSystemsBase.MatrixPencils.spzeros`](https://andreasvarga.github.io/MatrixPencils.jl/dev/sklfapps.html#MatrixPencils.spzeros) for more details. If `extra = Val(true)`, the function returns `z, iz, KRInfo` where `z` are the transmission zeros, information on the multiplicities of infinite zeros in `iz` and information on the Kronecker-structure in the KRInfo object. The number of infinite zeros is the sum of the components of iz.

To compute zeros of a system with non-BLAS floats, such as `BigFloat`, install and load the package `GenericSchur.jl` before calling `tzeros`.
"""
function tzeros(sys::TransferFunction{<:TimeEvolution,SisoRational{T}}; kwargs...) where {T}
    if issiso(sys)
        return tzeros(sys.matrix[1,1])
    else
        n,d = numpoly(sys), denpoly(sys)
        filter!(isfinite, MatrixPencils.rmzeros(n, d; kwargs...)[1]) # This uses rm2ls -> rm2lspm internally
    end
end

function tzeros(sys::TransferFunction{<:TimeEvolution,SisoZpk{T,TR}}) where {T, TR}
    if issiso(sys)
        return tzeros(sys.matrix[1,1])
    else
        return tzeros(ss(sys))  # Convert to ss since rmzeros does this anyways, so no reason to pass through tf{SisoRational}
    end
end

tzeros(sys::AbstractStateSpace; kwargs...) = tzeros(sys.A, sys.B, sys.C, sys.D; kwargs...)
# Make sure everything is BlasFloat
function tzeros(A::AbstractMatrix, B::AbstractMatrix, C::AbstractMatrix, D::AbstractMatrix)
    T = promote_type(eltype(A), eltype(B), eltype(C), eltype(D))
    A2, B2, C2, D2, _ = promote(A,B,C,D, fill(zero(T)/one(T),0,0)) # If Int, we get Float64
    tzeros(A2, B2, C2, D2)
end

function tzeros(A::AbstractMatrix{T}, B::AbstractMatrix{T}, C::AbstractMatrix{T}, D::AbstractMatrix{T}; extra::Val{E} = Val{false}(), balance=true, kwargs...) where {T <: BlasFloat, E}
    if balance
        A, B, C = balance_statespace(A, B, C)
    end

    (z, iz, KRInfo) = MatrixPencils.spzeros(A, I, B, C, D; kwargs...)
    if E
        return (z, iz, KRInfo)
    else
        return filter(isfinite, z)
    end
end

function tzeros(A::AbstractMatrix{T}, B::AbstractMatrix{T}, C::AbstractMatrix{T}, D::AbstractMatrix{T}; extra::Val{E} = Val{false}(), kwargs...) where {T <: Union{AbstractFloat,Complex{<:AbstractFloat}}, E}
    isempty(A) && return complex(T)[]

    # Balance the system
    A, B, C = balance_statespace(A, B, C; verbose=false)

    # Compute a good tolerance
    meps = 10*eps(real(T))*norm([A B; C D])

    # Step 1:
    A_r, B_r, C_r, D_r = reduce_sys(A, B, C, D, meps)

    # Step 2: (conjugate transpose should be avoided since single complex zeros get conjugated)
    A_rc, B_rc, C_rc, D_rc = reduce_sys(copy(transpose(A_r)), copy(transpose(C_r)), copy(transpose(B_r)), copy(transpose(D_r)), meps)

    # Step 3:
    # Compress cols of [C D] to [0 Df]
    mat = [C_rc D_rc]
    Wr = qr(mat').Q * I
    W = reverse(Wr, dims=2)
    mat = mat*W
    if fastrank(mat', meps) > 0
        nf = size(A_rc, 1)
        m = size(D_rc, 2)
        Af = ([A_rc B_rc] * W)[1:nf, 1:nf]
        Bf = ([Matrix{T}(I, nf, nf) zeros(nf, m)] * W)[1:nf, 1:nf]
        Z = schur(complex.(Af), complex.(Bf)) # Schur instead of eigvals to handle BigFloat
        zs = Z.values
        ControlSystemsBase._fix_conjugate_pairs!(zs) # Generalized eigvals does not return exact conj. pairs
    else
        zs = complex(T)[]
    end
    return zs
end

"""
    reduce_sys(A::AbstractMatrix, B::AbstractMatrix, C::AbstractMatrix, D::AbstractMatrix, meps::AbstractFloat)
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
        V = reverse(qr(Ctilde').Q * I, dims=2)
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
        Vm = [V zeros(T, n, m); zeros(T, m, n) Matrix{T}(I, m, m)] # I(m) is not used for type stability reasons (as of julia v1.7)
        if sigma > 0
            M = [A B; Cbar Dbar]
            Vs = [copy(V') zeros(T, n, sigma) ; zeros(T, sigma, n) Matrix{T}(I, sigma, sigma)]
        else
            M = [A B]
            Vs = copy(V')
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
    relative_gain_array(G, w::AbstractVector)
    relative_gain_array(G, w::Number)

Calculate the relative gain array of `G` at frequencies `w`. 
G(iω) .* pinv(tranpose(G(iω)))

The RGA can be used to find input-output pairings for MIMO control using individually tuned loops. Pair the inputs and outputs such that the RGA(ωc) at the crossover frequency becomes as close to diagonal as possible. Avoid pairings such that RGA(0) contains negative diagonal elements. 

- The sum of the absolute values of the entries in the RGA is a good measure of the "true condition number" of G, the best condition number that can be achieved by input/output scaling of `G`, -Glad, Ljung.
- The RGA is invariant to input/output scaling of `G`.
- If the RGA contains large entries, the system may be sensitive to model errors, -Skogestad, "Multivariable Feedback Control: Analysis and Design":
    - Uncertainty in the input channels (diagonal input uncertainty). Plants with
    large RGA-elements around the crossover frequency are fundamentally
    difficult to control because of sensitivity to input uncertainty (e.g. caused
    by uncertain or neglected actuator dynamics). In particular, decouplers or
    other inverse-based controllers should not be used for plants with large RGAeleme
    - Element uncertainty. Large RGA-elements imply sensitivity to element-by-element uncertainty.
    However, this kind of uncertainty may not occur in practice due to physical couplings
    between the transfer function elements. Therefore, diagonal input uncertainty
    (which is always present) is usually of more concern for plants with large RGA elements.

The relative gain array is computed using the The unit-consistent (UC) generalized inverse
Reference: "On the Relative Gain Array (RGA) with Singular and Rectangular Matrices"
Jeffrey Uhlmann
https://arxiv.org/pdf/1805.10312.pdf
"""
function relative_gain_array(G, w::AbstractVector)
    mapslices(relative_gain_array, freqresp(G, w), dims=(1,2))
end

relative_gain_array(G, w::Number) = relative_gain_array(freqresp(G, w))

"""
    relative_gain_array(A::AbstractMatrix; tol = 1.0e-15)

Reference: "On the Relative Gain Array (RGA) with Singular and Rectangular Matrices"
Jeffrey Uhlmann
https://arxiv.org/pdf/1805.10312.pdf
"""
function relative_gain_array(A::AbstractMatrix; tol = 1e-15)
    m, n = size(A)
    if m == n && LinearAlgebra.det(A) != 0
        return A .* inv(transpose(A))
    end
    L = zeros(m, n)
    M = ones(m, n)
    S = sign.(A)
    AA = abs.(A)
    idx = findall(vec(AA) .> 0.0)
    L[idx] = log.(AA[idx])
    idx = setdiff(1 : length(AA), idx)
    L[idx] .= 0
    M[idx] .= 0
    r = sum(M, dims=2)[:]
    c = sum(M, dims=1)[:]
    u = zeros(m, 1)
    v = zeros(1, n)
    dx = 2*tol
    while dx > tol
        idx = c .> 0
        p = sum(L[:, idx], dims=1) ./ c[idx]'
        L[:, idx] .-= p .* M[:, idx]
        v[idx] .-= p'
        dx = sum(abs, p)/length(p)
        idx = r .> 0
        p = sum(L[idx, :], dims=2) ./ r[idx]
        L[idx, :] .-= p .* M[idx, :]
        u[idx] .-= p
        dx += sum(abs, p)/length(p)
    end
    dl = exp.(u)
    dr = exp.(v)
    S = S .* exp.(L)
    A .* transpose(pinv(S) .* (dl * dr)')
end

"""
    wgm, gm, wpm, pm = margin(sys::LTISystem, w::Vector; full=false, allMargins=false, adjust_phase_start=true)

returns frequencies for gain margins, gain margins, frequencies for phase margins, phase margins

- If `!allMargins`, return only the smallest margin
- If `full` return also `fullPhase`
- `adjust_phase_start`: If true, the phase will be adjusted so that it starts at -90*intexcess degrees, where `intexcess` is the integrator excess of the system.

See also [`delaymargin`](@ref) and [`RobustAndOptimalControl.diskmargin`](https://juliacontrol.github.io/RobustAndOptimalControl.jl/dev/api/#RobustAndOptimalControl.diskmargin)
"""
function margin(sys::LTISystem, w::AbstractVector{<:Real}; full=false, allMargins=false, adjust_phase_start=true)
    ny, nu = size(sys)

    T = float(numeric_type(sys))
    if allMargins
        wgm         = Array{Array{T,1}}(undef, ny,nu)
        gm          = Array{Array{T,1}}(undef, ny,nu)
        wpm         = Array{Array{T,1}}(undef, ny,nu)
        pm          = Array{Array{T,1}}(undef, ny,nu)
        fullPhase   = Array{Array{T,1}}(undef, ny,nu)
    else
        wgm         = Array{T,2}(undef, ny, nu)
        gm          = Array{T,2}(undef, ny, nu)
        wpm         = Array{T,2}(undef, ny, nu)
        pm          = Array{T,2}(undef, ny, nu)
        fullPhase   = Array{T,2}(undef, ny, nu)
    end
    for j=1:nu
        for i=1:ny
            wgm[i,j], gm[i,j], wpm[i,j], pm[i,j], fullPhase[i,j] = sisomargin(sys[i,j], w; full=true, allMargins, adjust_phase_start)
        end
    end
    if full
        (; wgm, gm, wpm, pm, fullPhase)
    else
        (; wgm, gm, wpm, pm)
    end
end

"""
    ωgm, gm, ωpm, pm = sisomargin(sys::LTISystem, w::Vector; full=false, allMargins=false, adjust_phase_start=true))

Return frequencies for gain margins, gain margins, frequencies for phase margins, phase margins. If `allMargins=false`, only the smallest margins are returned.
"""
function sisomargin(sys::LTISystem, w::AbstractVector{<:Real}; full=false, allMargins=false, adjust_phase_start=true)
    ny, nu = size(sys)
    if ny !=1 || nu != 1
        error("System must be SISO, use `margin` instead")
    end
    mag, phase, w = bode(sys, w)
    wgm, = _allPhaseCrossings(w, phase)
    gm = similar(wgm)
    for i = eachindex(wgm)
        gm[i] = 1 ./ abs(freqresp(sys,wgm[i])[1])
    end
    wpm, fi = _allGainCrossings(w, mag)
    pm = similar(wpm)
    for i = eachindex(wpm)
        # We have to access the actual phase value from the `phase` array to get unwrapped phase. This value is not fully accurate since it is computed at a grid point, so we compute the more accurate phase at the interpolated frequency. This accurate value is not unwrapped, so we add an integer multiple of 360 to get the closest unwrapped phase.
        φ_nom = rad2deg(angle(freqresp(sys,wpm[i])[1]))
        φ_rounded = phase[clamp(round(Int, fi[i]), 1, length(phase))] # fi is interpolated, so we round to the closest integer
        φ_int = φ_nom - 360 * round( (φ_nom - φ_rounded) / 360 )

        # Now compute phase margin relative to -180:
        pm[i] = 180 + φ_int
    end
    if !allMargins #Only output the smallest margins
        gm, idx = findmin([gm;Inf])
        wgm = [wgm;NaN][idx]
        pm, idx = findmin([pm;Inf])
        fi = [fi;NaN][idx]
        wpm = [wpm;NaN][idx]
        if full
            if !isnan(fi) #fi may be NaN, fullPhase is a scalar
                fullPhase = interpolate(fi, phase)
            else
                fullPhase = NaN
            end
        end
    else
        isempty(gm) && (gm = [Inf])
        if full #We know that all values are defined and fullPhase is a vector
            fullPhase = interpolate(fi, phase)
        end
    end
    if adjust_phase_start && isrational(sys)
        intexcess = integrator_excess(sys)
        if intexcess != 0
            # Snap phase so that it starts at -90*intexcess
            nineties = round(Int, phase[1] / 90)
            adjust = ((90*(-intexcess-nineties)) ÷ 360) * 360
            pm = pm .+ adjust
            phase .+= adjust
            fullPhase = fullPhase .+ adjust
        end
    end
    if full
        (; wgm, gm, wpm, pm, fullPhase, phasedata = phase[:]) # phasedata is used by marginplot
    else
        (; wgm, gm, wpm, pm)
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
            push!(wcross, w[i])
            push!(tcross, i)
        elseif n[i] != n[i+1]
            #Interpolate to approximate crossing
            t = res[i]/(res[i]-res[i+1])
            push!(tcross, i+t)
            wt = w[i]+t*(w[i+1]-w[i])
            push!(wcross, wt)
        end
    end
    if res[end] == 0 #Special case if multiple points
        push!(wcross, w[end])
        push!(tcross, length(w))
    end
    wcross, tcross
end

"""
    dₘ = delaymargin(G::LTISystem)

Return the delay margin, dₘ. For discrete-time systems, the delay margin is normalized by the sample time, i.e., the value represents the margin in number of sample times. 
Only supports SISO systems.
"""
function delaymargin(G::LTISystem)
    # Phase margin in radians divided by cross-over frequency in rad/s.
    issiso(G) || error("delaymargin only supports SISO systems")
    m     = margin(G,allMargins=true)
    isempty(m[4][1]) && return Inf
    ϕₘ, i = findmin(m[4][1])
    ϕₘ   *= π/180
    ωϕₘ   = m[3][1][i]
    dₘ    = ϕₘ/ωϕₘ
    if isdiscrete(G)
        dₘ /= G.Ts # Give delay margin in number of sample times, as matlab does
    end
    dₘ
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
- `PS = (1+PC)\\P` Load disturbance to measurement signal
- `CS = (1+PC)\\C` Measurement noise to control signal
- `T = PC/(1+PC)` Complementary sensitivity function

If `minimal=true`, [`minreal`](@ref) will be applied to all transfer functions.
"""
function gangoffour(P::LTISystem, C::LTISystem; minimal=true)
    minfun = minimal ? robust_minreal : identity
    S = feedback(I(noutputs(P)), P*C)    |> minfun
    PS = feedback(P, C)     |> minfun
    CS = feedback(C, P)     |> minfun
    T = feedback(P*C, I(noutputs(P)))    |> minfun
    return S, PS, CS, T
end

"""
    S, PS, CS, T, RY, RU, RE = gangofseven(P,C,F)

Given transfer functions describing the Plant `P`, the controller `C` and a feed forward block `F`,
computes the four transfer functions in the Gang-of-Four and the transferfunctions corresponding to the feed forward.

- `S = 1/(1+PC)` Sensitivity function
- `PS = P/(1+PC)`
- `CS = C/(1+PC)`
- `T = PC/(1+PC)` Complementary sensitivity function
- `RY = PCF/(1+PC)`
- `RU = CF/(1+P*C)`
- `RE = F/(1+P*C)`
"""
function gangofseven(P::LTISystem, C::LTISystem, F::LTISystem)
    S, PS, CS, T = gangoffour(P,C)
    RY = T*F
    RU = CS*F
    RE = S*F
    return S, PS, CS, T, RY, RU, RE
end
