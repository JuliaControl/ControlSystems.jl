function freqresp!(R::Array{T,3}, sys::DelayLtiSystem, ω::AbstractVector{W}) where {T, W <: Real}
    ny = noutputs(sys)
    nu = ninputs(sys)
    @boundscheck size(R) == (ny,nu,length(ω))
    P_fr = freqresp(sys.P.P, ω)

    cache = cis.(ω[1].*sys.Tau)

    @views for ω_idx = eachindex(ω)
        P11_fr = P_fr[1:ny, 1:nu, ω_idx]
        P12_fr = P_fr[1:ny, nu+1:end, ω_idx]
        P21_fr = P_fr[ny+1:end, 1:nu, ω_idx]
        P22_fr = P_fr[ny+1:end, nu+1:end, ω_idx]
        @. cache = cis(ω[ω_idx]*sys.Tau)
        delay_matrix_inv_fr = Diagonal(cache) # Frequency response of the diagonal matrix with delays
        # Inverse of the delay matrix, so there should not be any minus signs in the exponents

        R[:,:,ω_idx] .= P11_fr .+ P12_fr/(delay_matrix_inv_fr - P22_fr)*P21_fr # The matrix is invertible (?!)
    end
    return R
end

function evalfr(sys::DelayLtiSystem, s)
    (ny, nu) = size(sys)

    P_fr = evalfr(sys.P.P, s)

    P11_fr = P_fr[1:ny, 1:nu]
    P12_fr = P_fr[1:ny, nu+1:end]
    P21_fr = P_fr[ny+1:end, 1:nu]
    P22_fr = P_fr[ny+1:end, nu+1:end]

    delay_matrix_inv_fr = Diagonal(exp.(s*sys.Tau))

    return P11_fr + P12_fr/(delay_matrix_inv_fr - P22_fr)*P21_fr
end


"""
    delayd_ss(τ, Ts)

Discrete-time statespace realization of a delay ``τ`` sampled with period ``T_s``,
i.e. of ``z^{-N}`` where ``N = τ / T_s.``

``τ`` must be a multiple of ``T_s``. See [`thiran`](@ref) for approximate discretization of fractional delays.
"""
function delayd_ss(τ::Number, Ts::Number)
    n = round(Int, τ / Ts)
    if !(τ - n*Ts ≈ 0)
        error("The delay τ must be a multiple of the sample time Ts, use the function `thiran` to approximately discretize fractional delays.")
    end
    ss(diagm(1 => ones(n-1)), [zeros(n-1,1); 1], [1 zeros(1,n-1)], 0, Ts)
end

"""
    c2d(G::DelayLtiSystem, Ts, method=:zoh)
"""
function c2d(G::DelayLtiSystem, Ts::Real, method=:zoh)
    if !(method === :zoh)
        error("c2d for DelayLtiSystems only supports zero-order hold")
    end
    X = append(delayd_ss(τ, Ts) for τ in G.Tau)
    Pd = c2d(G.P.P, Ts)
    return lft(Pd, X)
end


# We have to default to something, look at the sys.P.P and delays
function _bounds_and_features(sys::DelayLtiSystem, plot::Val)
    ws, pz =  _bounds_and_features(sys.P.P, plot)
    logtau = log10.(abs.(sys.Tau))
    logtau = filter(x->x>4, logtau) # Ignore low frequency
    if isempty(logtau)
        return ws, pz
    end
    extreme = extrema(logtau)
    return [min(ws[1], floor(extreme[1]-0.2)), max(ws[2], ceil(extreme[2]+0.2))], pz
end

# Again we have to do something for default vectors, more or less a copy from timeresp.jl
function _default_dt(sys::DelayLtiSystem)
    if !isstable(sys.P.P)
        return 0.05   # Something small
    else
        ps = poles(sys.P.P)
        r = minimum([abs.(real.(ps));0]) # Find the fastest pole of sys.P.P
        r = min(r, minimum([sys.Tau;0])) # Find the fastest delay
        if r == 0.0
            r = 1.0
        end
        return 0.07/r
    end
end


# Used for pade approximation
"""
    p2 = _linscale(p::Polynomial, a)

Given a polynomial `p` and a number `a, returns the polynomials `p2` such that
`p2(s) == p(a*s)`.
"""
function _linscale(p::Polynomial, a)
    # This function should perhaps be implemented in Polynomials.jl
    coeffs_scaled = similar(p.coeffs, promote_type(eltype(p), typeof(a)))
    a_pow = 1
    coeffs_scaled[1] = p.coeffs[1]
    for k=2:length(p.coeffs)
        a_pow *= a
        coeffs_scaled[k] = p.coeffs[k]*a_pow
    end
    return Polynomial(coeffs_scaled)
end

# Coefficients for Padé approximations
# PADE_Q_COEFFS = [Polynomial([binomial(N,i)*prod(N+1:2*N-i) for i=0:N]) for N=1:10]
const PADE_Q_COEFFS =  [[2, 1],
 [12, 6, 1],
 [120, 60, 12, 1],
 [1680, 840, 180, 20, 1],
 [30240, 15120, 3360, 420, 30, 1],
 [665280, 332640, 75600, 10080, 840, 42, 1],
 [17297280, 8648640, 1995840, 277200, 25200, 1512, 56, 1],
 [518918400, 259459200, 60540480, 8648640, 831600, 55440, 2520, 72, 1],
 [17643225600, 8821612800, 2075673600, 302702400, 30270240, 2162160, 110880, 3960, 90, 1],
 [670442572800, 335221286400, 79394515200, 11762150400, 1210809600, 90810720, 5045040, 205920, 5940, 110, 1]]

"""
    pade(τ::Real, N::Int)

Compute the `N`th order Padé approximation of a time-delay of length `τ`.

See also [`thiran`](@ref) for discretization of delays.
"""
function pade(τ::Real, N::Int)
    if !(1 <= N <= 10); error("Order of Padé approximation must be between 1 and 10. Got $N."); end

    Q = Polynomials.Polynomial(PADE_Q_COEFFS[N])

    return tf(_linscale(Q, -τ), _linscale(Q, τ)) # return Q(-τs)/Q(τs)
end


"""
    pade(G::DelayLtiSystem, N)

Approximate all time-delays in `G` by Padé approximations of degree `N`.
"""
function pade(G::DelayLtiSystem, N)
    ny, nu = size(G)
    nTau = length(G.Tau)
    X = append(ss(pade(τ,N)) for τ in G.Tau) # Perhaps append should be renamed blockdiag
    return lft(G.P.P, X)
end

"""
    thiran(τ::Real, Ts)

Discretize a potentially fractional delay ``τ`` as a Thiran all-pass filter with sample time `Ts`. 

The Thiran all-pass filter gives an a maximally flat group delay.

If ``τ`` is an integer multiple of ``Ts``, the Thiran all-pass filter reduces to ``z^{-τ/Ts}``.

Ref: T. I. Laakso, V. Valimaki, M. Karjalainen and U. K. Laine, "Splitting the unit delay [FIR/all pass filters design]," in IEEE Signal Processing Magazine, vol. 13, no. 1, 1996.
"""
function thiran(τ::Real, Ts)
    D = τ/Ts
    N = ceil(Int, D)
    a = ones(N+1)
    for k = 1:N
        P = prod((D-N+n) / (D-N+n+k) for n in 0:N)
        a[k+1] = (-1)^k * binomial(N, k) * P # Eq 86 in reference
    end
    tf(reverse(a), a, Ts)
end