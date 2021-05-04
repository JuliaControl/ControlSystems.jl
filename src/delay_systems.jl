function freqresp(sys::DelayLtiSystem, ω::AbstractVector{T}) where {T <: Real}
    ny = noutputs(sys)
    nu = ninputs(sys)

    P_fr = freqresp(sys.P.P, ω);

    G_fr = zeros(eltype(P_fr), length(ω), ny, nu)

    for ω_idx=1:length(ω)
        P11_fr = P_fr[ω_idx, 1:ny, 1:nu]
        P12_fr = P_fr[ω_idx, 1:ny, nu+1:end]
        P21_fr = P_fr[ω_idx, ny+1:end, 1:nu]
        P22_fr = P_fr[ω_idx, ny+1:end, nu+1:end]

        delay_matrix_inv_fr = Diagonal(exp.(im*ω[ω_idx]*sys.Tau)) # Frequency response of the diagonal matrix with delays
        # Inverse of the delay matrix, so there should not be any minus signs in the exponents

        G_fr[ω_idx,:,:] .= P11_fr + P12_fr/(delay_matrix_inv_fr - P22_fr)*P21_fr # The matrix is invertible (?!)
    end

    return G_fr
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
Discrete-time statespace realization of a delay τ sampled with period Ts,
i.e. of z^-N where N = τ / Ts.

τ must be a multiple of Ts.
"""
function delayd_ss(τ::Number, Ts::Number)
    n = Int(round(τ / Ts))
    if !(τ - n*Ts ≈ 0)
        error("The delay τ must be a multiple of the sample time Ts")
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
    X = append([delayd_ss(τ, Ts) for τ in G.Tau]...)
    Pd = c2d(G.P.P, Ts)
    return lft(Pd, X)
end


"""
    `y, t, x = lsim(sys::DelayLtiSystem, u, t::AbstractArray{<:Real}; x0=fill(0.0, nstates(sys)), alg=MethodOfSteps(Tsit5()), abstol=1e-6, reltol=1e-6, kwargs...)`

    Simulate system `sys`, over time `t`, using input signal `u`, with initial state `x0`, using method `alg` .

    Arguments:

    `t`: Has to be an `AbstractVector` with equidistant time samples (`t[i] - t[i-1]` constant)
    `u`: Function to determine control signal `ut` at a time `t`, on any of the following forms:
        Can be a constant `Number` or `Vector`, interpreted as `ut .= u` , or
        Function `ut .= u(t)`, or
        In-place function `u(ut, t)`. (Slightly more effienct)
    `alg, abstol, reltol` and `kwargs...`: are sent to `DelayDiffEq.solve`.

    Returns: times `t`, and `y` and `x` at those times.
"""
function lsim(sys::DelayLtiSystem{T,S}, u, t::AbstractArray{<:Real};
        x0=fill(zero(T), nstates(sys)),
        alg=DelayDiffEq.MethodOfSteps(Tsit5()),
        abstol=1e-6, reltol=1e-6,
        kwargs...) where {T,S}

    # Make u! in-place function of u
    u! = if isa(u, Number) || isa(u,AbstractVector) # Allow for u to be a constant number or vector
        (uout, t) -> uout .= u
    elseif DiffEqBase.isinplace(u, 2)               # If u is an inplace (more than 1 argument function)
        u
    else                                            # If u is a regular u(t) function
        (out, t) -> (out .= u(t))
    end

    _lsim(sys, u!, t, x0, alg; abstol=abstol, reltol=reltol, kwargs...)
end

# Generic parametrized dde used for simulating DelayLtiSystem
# We simulate the integral of the states
# u, du are from the notation of variables in DifferentialEquations.jl
# The state from the control system is x
# u!, uout is the control law and its output
function dde_param(du, u, h, p, t)
    A, B1, B2, C1, C2, D11, D12, D21, D22, y, Tau, u!, uout, hout, tmpy, tsave = p

    nx = size(A,1)
    nd = length(Tau)
    ny = size(C1,1)
    
    dx = view(du, 1:nx)
    dY = view(du, (nx+1):(nx+ny))
    dD = view(du, (nx+ny+1):(nx+ny+nd))
    x = view(u, 1:nx)

    # uout = u(t)
    u!(uout, t)

    # hout = d(t-Tau)
    for k=1:length(Tau)
        # Get the derivative of the history for d at indices corresponding to Tau
        hout[k] = h(p, t-Tau[k], Val{1}, idxs=(nx+ny+k))
    end

    #dx(t) .= A*x(t) + B1*u(t) +B2*d(t-tau)
    mul!(dx, A, x)
    mul!(dx, B1, uout, true, true)
    mul!(dx, B2, hout, true, true)

    # dY = y(t) = C1*x + D11*u(t) + D12*d(t-Tau)
    mul!(dY, C1, x)
    mul!(dY, D11, uout, true, true)
    mul!(dY, D12, hout, true, true)

    # dD = d(t) = C2*x + D21*u(t) + D22*d(t-Tau)
    mul!(dD, C2, x)
    mul!(dD, D21, uout, true, true)
    mul!(dD, D22, hout, true, true)

    # Save y value in tmpy to be used by dde_saver
    if t in tsave
        # The value of y at this time is given by the derivative
        tmpy .= dY
    end
    return
end

# Save x(t) and y(t) to output
function dde_saver(u,t,integrator)
    A, B1, B2, C1, C2, D11, D12, D21, D22, y, Tau, u!, uout, hout, tmpy, tsave = integrator.p
    nx = size(A,1)
    # y is already saved in tmpy
    u[1:nx], copy(tmpy)
end

function _lsim(sys::DelayLtiSystem{T,S}, Base.@nospecialize(u!), t::AbstractArray{<:Real}, x0::Vector{T}, alg; kwargs...) where {T,S}

    P = sys.P

    t0 = first(t)
    dt = t[2] - t[1]

    # Get all matrices to save on allocations
    A, B1, B2, C1, C2, D11, D12, D21, D22 = P.A, P.B1, P.B2, P.C1, P.C2, P.D11, P.D12, P.D21, P.D22
    Tau = sys.Tau
    nx = size(A,1)
    nd = length(Tau)
    ny = size(C1,1)
    nt = length(t)

    hout = fill(zero(T), nd)            # in place storage for delays
    uout = fill(zero(T), ninputs(sys))  # in place storage for u
    tmpy = similar(x0, ny)              # in place storage for output
    y = Matrix{T}(undef, ny, nt)        # Output matrix
    x = Matrix{T}(undef, nx, nt)        # State matrix

    p = (A, B1, B2, C1, C2, D11, D12, D21, D22, y, Tau, u!, uout, hout, tmpy, t)

    # This callback computes and stores the delay term
    sv = SavedValues(Float64, Tuple{Vector{Float64},Vector{Float64}})
    cb = SavingCallback(dde_saver, sv, saveat = t)

    # History function, only used for d
    # The true d(t) is stored as a derivative
    h!(p, t, deriv::Type{Val{1}}; idxs=0) = zero(T)

    # The states are x(t), Y(t), D(t), where Y, D are integrals of y(t), d(t)
    u0 = [x0;zeros(T,ny);zeros(T,nd)]
    prob = DDEProblem{true}(dde_param, u0, h!,
                (T(t[1]), T(t[end])),
                p,
                constant_lags=sort(Tau),# Not sure if sort needed
                neutral=true,           # We have derivatives on RHS (d(t)) 
                callback=cb)
    # Important to stop at t since we can not access derivatives in SavingCallback
    sol = DelayDiffEq.solve(prob, alg; tstops=t, saveat=t, kwargs...)

    # Retrive the saved values
    for k = 1:nt
        x[:,k] .= sv.saveval[k][1]
        y[:,k] .= sv.saveval[k][2]
    end

    return y', t, x'
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

# Againm we have to do something for default vectors, more or less a copy from timeresp.jl
function _default_dt(sys::DelayLtiSystem)
    if !isstable(sys.P.P)
        return 0.05   # Something small
    else
        ps = pole(sys.P.P)
        r = minimum([abs.(real.(ps));0]) # Find the fastest pole of sys.P.P
        r = min(r, minimum([sys.Tau;0])) # Find the fastest delay
        if r == 0.0
            r = 1.0
        end
        return 0.07/r
    end
end

function Base.step(sys::DelayLtiSystem{T}, t::AbstractVector; kwargs...) where T
    nu = ninputs(sys)
    if t[1] != 0
        throw(ArgumentError("First time point must be 0 in step"))
    end
    u = (out, t) -> (t < 0 ? out .= 0 : out .= 1)
    x0=fill(zero(T), nstates(sys))
    if nu == 1
        y, tout, x = lsim(sys, u, t; x0=x0, kwargs...)
    else
        x = Array{T}(undef, length(t), nstates(sys), nu)
        y = Array{T}(undef, length(t), noutputs(sys), nu)
        for i=1:nu
            y[:,:,i], tout, x[:,:,i] = lsim(sys[:,i], u, t; x0=x0, kwargs...)
        end
    end
    return y, tout, x
end


function impulse(sys::DelayLtiSystem{T}, t::AbstractVector; alg=MethodOfSteps(BS3()), kwargs...) where T
    nu = ninputs(sys)
    iszero(sys.P.D12) || @warn("Impulse with a direct term from input to delay vector leads to poor accuracy.")
    if t[1] != 0
        throw(ArgumentError("First time point must be 0 in impulse"))
    end
    u = (out, t) -> (out .= 0)
    if nu == 1
        y, tout, x = lsim(sys, u, t; alg=alg, x0=sys.P.B[:,1], kwargs...)
    else
        x = Array{T}(undef, length(t), nstates(sys), nu)
        y = Array{T}(undef, length(t), noutputs(sys), nu)
        for i=1:nu
            y[:,:,i], tout, x[:,:,i] = lsim(sys[:,i], u, t; alg=alg, x0=sys.P.B[:,i], kwargs...)
        end
    end
    return y, tout, x
end


# Used for pade approximation
"""
`p2 = _linscale(p::Polynomial, a)`

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

# Coefficeints for Padé approximations
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
    X = append([ss(pade(τ,N)) for τ in G.Tau]...) # Perhaps append should be renamed blockdiag
    return lft(G.P.P, X)
end
