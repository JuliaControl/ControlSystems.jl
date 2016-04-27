@doc """`pole(sys)`

Compute the poles of system `sys`.""" ->
pole(sys::StateSpace) = eig(sys.A)[1]
pole(sys::TransferFunction) = [map(pole, sys.matrix)...;]
pole(sys::SisoTf) = error("pole is not implemented for type $(typeof(sys))")

@doc """`dcgain(sys)`

Compute the gain of SISO system `sys`.""" ->
function dcgain(sys::StateSpace, zs::Vector=tzero(sys))
    !issiso(sys) && error("Gain only defined for siso systems")
    return ( sys.C*(-sys.A\sys.B) + sys.D)[1]
end

# Helper function for SISO systems
function _dcgain(sys::SisoTf, discrete::Bool = false)
  zerovec = numvec(sys)[end:-1:1]
  polevec = denvec(sys)[end:-1:1]
  K       = length(polevec)

  # Assume causal system. `numvec` and `denvec` functions do not return
  # vectors of the same size. That's why try to have the same sized arrays.
  zerovec = vcat( zerovec, zeros(K-length(zerovec)) )

  if discrete   # Discrete systems
    # Right now `isstable` does not work with sys::SisoTf inputs.
    # if !isstable(sys)
    #   return NaN
    # end

    for k = 1:K
      zerosum = sum(zerovec[k:end])
      polesum = sum(polevec[k:end])

      # 0/0 will result in L'Hopital
      if zerosum == polesum == zero(zerosum)
        zerovec[k+1:end] = zerovec[k+1:end].*(1:K-k)
        polevec[k+1:end] = polevec[k+1:end].*(1:K-k)
      else
        return zerosum/polesum
      end
    end
  else          # Continuous systems
    # Right now `isstable` does not work with sys::SisoTf inputs.
    # if !isstable(sys)
    #   return NaN
    # end

    for k = 1:K
      # 0/0 will result in L'Hopital
      if zerovec[k] == polevec[k] == zero(zerovec[1])
        zerovec[k+1:end] = zerovec[k+1:end].*(1:K-k)
        polevec[k+1:end] = polevec[k+1:end].*(1:K-k)
      else
        return zerovec[k]/polevec[k]
      end
    end
  end
end

"""
    dcgain(sys::TransferFunction)

Return the DC gain of system `sys`.

The function returns `NaN` for unstable systems, and ``G(s=0)`` and ``G(z=1)`` for stable
continuous- and discrete-time systems, respectively.

# Examples
```julia
julia> # Continuous-time systems

julia> c1 = tf([1],[1,0]); # 1/s

julia> c2 = tf([1,0],[1,1]); # s/(s+1)

julia> c3 = tf([1],[1,5]); # s/(s+5)

julia> c4 = tf([1],[1,-5]); # s/(s-5)

julia> c5 = tf([c1 c2;c3 c4]);

julia> dcgain(c5)
2x2 Array{Float64,2}:
 Inf      0.0
   0.2  NaN  

julia> # Discrete-time systems

julia> d1 = tf([1],[1,-1],0.5); # 1/(z-1)

julia> d2 = tf([1,-1],[1,0],0.5); # (z-1)/z

julia> d3 = tf([1,-.7],[1,-.3],0.5); # (z-.7)/(z-.3)

julia> d4 = tf([1,-.7],[1,-1.5],0.5); # (z-.7)/(z-1.5)

julia> d5 = tf([d1 d2;d3 d4]);

julia> dcgain(d5)
2x2 Array{Float64,2}:
 Inf           0.0
    0.428571  NaN  
```
"""
function dcgain(sys::TransferFunction)
  gains = zeros(size(sys.matrix))

  for idx in eachindex(sys.matrix)
    gains[idx] = _dcgain(sys.matrix[idx], sys.Ts > 0)
  end

  return gains
end

@doc """`markovparam(sys, n)`

Compute the `n`th markov parameter of state-space system `sys`. This is defined
as the following:

`h(0) = D`

`h(n) = C*A^(n-1)*B`""" ->
function markovparam(sys::StateSpace, n::Integer)
    n < 0 && error("n must be >= 0")
    return n == 0 ? sys.D : sys.C * sys.A^(n-1) * sys.B
end

@doc """`z, p, k = zpkdata(sys)`

Compute the zeros, poles, and gains of system `sys`.

### Returns
`z` : Matrix{Vector{Complex128}}, (ny x nu)

`p` : Matrix{Vector{Complex128}}, (ny x nu)

`k` : Matrix{Float64}, (ny x nu)""" ->
function zpkdata(sys::LTISystem)
    ny, nu = size(sys)
    zs = Array(Vector{Complex128}, ny, nu)
    ps = Array(Vector{Complex128}, ny, nu)
    ks = Array(Float64, ny, nu)
    for j = 1:nu
        for i = 1:ny
            zs[i, j], ps[i, j], ks[i, j] = _zpk_kern(sys, i, j)
        end
    end
    return zs, ps, ks
end

function _zpk_kern(sys::StateSpace, iy::Int, iu::Int)
    A, B, C = struct_ctrb_obsv(sys.A, sys.B[:, iu:iu], sys.C[iy:iy, :])
    D = sys.D[iy:iy, iu:iu]
    z = tzero(A, B, C, D)
    nx = size(A, 1)
    nz = length(z)
    k = nz == nx ? D[1] : (C*(A^(nx - nz - 1))*B)[1]
    return z, eigvals(A), k
end

function _zpk_kern(sys::TransferFunction, iy::Int, iu::Int)
    s = sys.matrix[iy, iu]
    return _zpk_kern(s)
end

function _zpk_kern(s::SisoRational)
  return roots(s.num), roots(s.den), s.num[1]/s.den[1]
end

function _zpk_kern(s::SisoZpk)
    return s.z, s.p, s.k
end

@doc """`Wn, zeta, ps = damp(sys)`

Compute the natural frequencies, `Wn`, and damping ratios, `zeta`, of the
poles, `ps`, of `sys`""" ->
function damp(sys::LTISystem)
    ps = pole(sys)
    if !iscontinuous(sys)
        Ts = sys.Ts == -1 ? 1 : sys.Ts
        ps = log(ps)/Ts
    end
    Wn = abs(ps)
    order = sortperm(Wn)
    Wn = Wn[order]
    ps = ps[order]
    ζ = -cos(angle(ps))
    return Wn, ζ, ps
end

@doc """`dampreport(sys)`

Display a report of the poles, damping ratio, natural frequency, and time
constant of the system `sys`""" ->
function dampreport(io::IO, sys::LTISystem)
    Wn, zeta, ps = damp(sys)
    t_const = 1./(Wn.*zeta)
    header =
    ("|     Pole      |   Damping     |   Frequency   | Time Constant |\n"*
    "|               |    Ratio      |   (rad/sec)   |     (sec)     |\n"*
    "+---------------+---------------+---------------+---------------+")
    println(io, header)
    for i=eachindex(ps)
        p, z, w, t = ps[i], zeta[i], Wn[i], t_const[i]
        @printf(io, "|  %-13.3e|  %-13.3e|  %-13.3e|  %-13.3e|\n", p, z, w, t)
    end
end
dampreport(sys::LTISystem) = dampreport(STDOUT, sys)


@doc """`tzero(sys)`

Compute the invariant zeros of the system `sys`. If `sys` is a minimal
realization, these are also the transmission zeros.""" ->
function tzero(sys::TransferFunction)
    if issiso(sys)
        return tzero(sys.matrix[1,1])
    else
        return tzero(ss(sys))
    end
end

# Implements the algorithm described in:
# Emami-Naeini, A. and P. Van Dooren, "Computation of Zeros of Linear
# Multivariable Systems," Automatica, 18 (1982), pp. 415–430.
#
# Note that this returns either Vector{Complex64} or Vector{Float64}
tzero(sys::StateSpace) = tzero(sys.A, sys.B, sys.C, sys.D)
function tzero(A::Matrix{Float64}, B::Matrix{Float64}, C::Matrix{Float64},
        D::Matrix{Float64})
    # Balance the system
    A, B, C = balance_statespace(A, B, C)

    # Compute a good tolerance
    meps = 10*eps()*norm([A B; C D])
    A, B, C, D = reduce_sys(A, B, C, D, meps)
    A, B, C, D = reduce_sys(A', C', B', D', meps)
    if isempty(A)   return Float64[]    end

    # Compress cols of [C D] to [0 Df]
    mat = [C D]
    # To ensure type-stability, we have to annote the type here, as qrfact
    # returns many different types.
    W = full(qrfact(mat')[:Q], thin=false)::Matrix{Float64}
    W = flipdim(W,2)
    mat = mat*W
    if fastrank(mat', meps) > 0
        nf = size(A, 1)
        m = size(D, 2)
        Af = ([A B] * W)[1:nf, 1:nf]
        Bf = ([eye(nf) zeros(nf, m)] * W)[1:nf, 1:nf]
        zs = eig(Af, Bf)[1]
    else
        zs = Float64[]
    end
    return zs
end

"""
Implements REDUCE in the Emami-Naeini & Van Dooren paper. Returns transformed
A, B, C, D matrices. These are empty if there are no zeros.
"""
function reduce_sys(A::Matrix{Float64}, B::Matrix{Float64}, C::Matrix{Float64}, D::Matrix{Float64}, meps::Float64)
    Cbar, Dbar = C, D
    if isempty(A)
        return A, B, C, D
    end
    while true
        # Compress rows of D
        U = full(qrfact(D)[:Q], thin=false)::Matrix{Float64}
        D = U'*D
        C = U'*C
        sigma = fastrank(D, meps)
        Cbar = C[1:sigma, :]
        Dbar = D[1:sigma, :]
        Ctilde = C[(1 + sigma):end, :]
        if sigma == size(D, 1)
            break
        end

        # Compress columns of Ctilde
        V = full(qrfact(Ctilde')[:Q], thin=false)::Matrix{Float64}
        V = flipdim(V,2)
        Sj = Ctilde*V
        rho = fastrank(Sj', meps)
        nu = size(Sj, 2) - rho

        if rho == 0
            break
        elseif nu == 0
            # System has no zeros, return empty matrices
            A = B = Cbar = Dbar = Float64[]
            break
        end
        # Update System
        n, m = size(B)
        Vm = [V zeros(n, m); zeros(m, n) eye(m)]
        if sigma > 0
            M = [A B; Cbar Dbar]
            Vs = [V' zeros(n, sigma) ; zeros(sigma, n) eye(sigma)]
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
function fastrank(A::Matrix{Float64}, meps::Float64)
    n, m = size(A, 1, 2)
    if n*m == 0     return 0    end
    norms = Array(Float64, n)
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

"""
function margin{S<:Real}(sys::LTISystem, w::AbstractVector{S}; full=false, allMargins=false)
    ny, nu = size(sys)
    vals = (:wgm, :gm, :wpm, :pm, :fullPhase)
    if allMargins
        for val in vals
            eval(:($val = Array{Array{Float64,1}}($ny,$nu)))
        end
    else
        for val in vals
            eval(:($val = Array{Float64,2}($ny,$nu)))
        end
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

returns frequencies for gain margins, gain margins, frequencies for phase margins, phase margins
"""
function sisomargin{S<:Real}(sys::LTISystem, w::AbstractVector{S}; full=false, allMargins=false)
    ny, nu = size(sys)
    if ny !=1 || nu != 1
        error("System must be SISO, use `margin` instead")
    end
    mag, phase, w = bode(sys, w)
    wgm, = _allPhaseCrossings(w, phase)
    gm = similar(wgm)
    for i = eachindex(wgm)
        gm[i] = 1./abs(freqresp(sys,[wgm[i]])[1][1])
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
        pm, idx = findmin([abs(pm);Inf])
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
margin(system, _default_freq_vector(system, :bode); kwargs...)
#margin(sys::LTISystem, args...) = margin(LTISystem[sys], args...)

# Interpolate the values in "list" given the floating point "index" fi
function interpolate(fi, list)
    fif = floor(Integer, fi)
    fic = ceil(Integer, fi)
    list[fif]+mod(fi,1).*(list[fic]-list[fif])
end

function _allGainCrossings(w, mag)
    _findCrossings(w,mag.>1,mag-1)
end

function _allPhaseCrossings(w, phase)
    #Calculate numer of times real axis is crossed on negative side
    n =  Array{Float64,1}(length(w)) #Nbr of crossed
    ph = Array{Float64,1}(length(w)) #Residual
    for i = eachindex(w) #Found no easier way to do this
        n[i], ph[i] = fldmod(phase[i]+180,360)#+180
    end
    _findCrossings(w, n, ph)
end

function _findCrossings(w, n, res)
    wcross = Array{Float64,1}()
    tcross = Array{Float64,1}()
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

@doc """`dₘ = delaymargin(G::LTISystem)`

Only supports SISO systems""" ->
function delaymargin(G::LTISystem)
    # Phase margin in radians divided by cross-over frequency in rad/s.
    if G.nu + G.ny > 2
        error("delaymargin only supports SISO systems")
    end
    m   = margin(G,allMargins=true)
    ϕₘ,i= findmin(m[4])
    ϕₘ *= π/180
    ωϕₘ = m[3][i]
    dₘ  = ϕₘ/ωϕₘ
    if G.Ts > 0
        dₘ /= G.Ts # Give delay margin in number of sample times, as matlab does
    end
    dₘ
end

@doc """`S,D,N,T = gangoffour(P,C)`, `gangoffour(P::AbstractVector,C::AbstractVector)`

Given a transfer function describing the Plant `P` and a transferfunction describing the controller `C`, computes the four transfer functions in the Gang-of-Four.

`S = 1/(1+PC)` Sensitivity function

`D = P/(1+PC)`

`N = C/(1+PC)`

`T = PC/(1+PC)` Complementary sensitivity function

Only supports SISO systems""" ->
function gangoffour(P::TransferFunction,C::TransferFunction)
    if P.nu + P.ny + C.nu + C.ny > 4
        error("gangoffour only supports SISO systems")
    end
    S = minreal(1/(1+P*C))
    D = minreal(P*S)
    N = minreal(C*S)
    T = minreal(P*N)
    return S, D, N, T
end

function gangoffour(P::AbstractVector, C::AbstractVector)
    if P[1].nu + P[1].ny + C[1].nu + C[1].ny > 4
        error("gangoffour only supports SISO systems")
    end
    length(P) == length(C) || error("P has to be the same length as C")
    n = length(P)
    S = [minreal(1/(1+P[i]*C[i])) for i in 1:n]
    D = [minreal(P[i]*S[i]) for i in 1:n]
    N = [minreal(C[i]*S[i]) for i in 1:n]
    T = [minreal(P[i]*N[i]) for i in 1:n]
    return S, D, N, T
end

function gangoffour(P::TransferFunction, C::AbstractVector)
    gangoffour(fill(P,length(C)), C)
end

function gangoffour(P::AbstractVector, C::TransferFunction)
    gangoffour(P, fill(C,length(P)))
end

@doc """`S, D, N, T, RY, RU, RE = gangofseven(P,C,F)`

Given transfer functions describing the Plant `P`, the controller `C` and a feed forward block `F`,
computes the four transfer functions in the Gang-of-Four and the transferfunctions corresponding to the feed forward.

`S = 1/(1+PC)` Sensitivity function

`D = P/(1+PC)`

`N = C/(1+PC)`

`T = PC/(1+PC)` Complementary sensitivity function

`RY = PCF/(1+PC)`

`RU = CF/(1+P*C)`

`RE = F/(1+P*C)`

Only supports SISO systems""" ->
function gangofseven(P::TransferFunction,C::TransferFunction,F::TransferFunction)
    if P.nu + P.ny + C.nu + C.ny + F.nu + F.ny > 6
        error("gof only supports SISO systems")
    end
    S, D, N, T = gangoffour(P,C)
    RY = T*F
    RU = N*F
    RE = S*F
    return S, D, N, T, RY, RU, RE
end
