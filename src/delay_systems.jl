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
function delayd_ss(τ::Number, h::Number)
    n = Int(round(τ / h))
    if !(τ - n*h ≈ 0)
        error("The delay τ must be a multiple of the sample time Ts")
    end
    ss(diagm(1 => ones(n-1)), [zeros(n-1,1); 1], [1 zeros(1,n-1)], 0, h)
end

"""
    c2d(G::DelayLtiSystem, Ts, method=:zoh)
"""
function c2d(G::DelayLtiSystem, h::Real, method=:zoh)
    if !(method === :zoh)
        error("c2d for DelayLtiSystems only supports zero-order hold")
    end
    X = append([delayd_ss(τ, h) for τ in G.Tau]...)
    Pd = c2d(G.P.P, h)[1]
    return lft(Pd, X)
end


"""
    `y, t, x = lsim(sys::DelayLtiSystem, u, t::AbstractArray{<:Real}; x0=fill(0.0, nstates(sys)), alg=MethodOfSteps(Tsit5()), kwargs...)`

    Simulate system `sys`, over time `t`, using input signal `u`, with initial state `x0`, using method `alg` .

    Arguments:

    `t`: Has to be an `AbstractVector` with equidistant time samples (`t[i] - t[i-1]` constant)
    `u`: Function to determine control signal `ut` at a time `t`, on any of the following forms:
        Can be a constant `Number` or `Vector`, interpreted as `ut .= u` , or
        Function `ut .= u(t)`, or
        In-place function `u(ut, t)`. (Slightly more effienct)
    `kwargs...`: These are sent to `solve` from DelayDiffEq.

    Returns: times `t`, and `y` and `x` at those times.
"""
function lsim(sys::DelayLtiSystem{T,S}, u, t::AbstractArray{<:Real}; x0=fill(zero(T), nstates(sys)), alg=MethodOfSteps(Tsit5()), method= :simple, kwargs...) where {T,S}
    # Make u! in-place function of u
    u! = if isa(u, Number) || isa(u,AbstractVector) # Allow for u to be a constant number or vector
        (uout, t) -> uout .= u
    elseif DiffEqBase.isinplace(u, 2)               # If u is an inplace (more than 1 argument function)
        u
    else                                            # If u is a regular u(t) function
        (out, t) -> (out .= u(t))
    end

    # TODO MethodOfSteps?
    _lsim(sys, u!, t, x0, Tsit5(); kwargs...)
end

function dde_param(dx, x, p, t)
    A, B1, B2, C1, C2, D11, D12, D21, D22, y, Tau, u!, h, uout, hout, tmpx, tmpy1, tmpy2, tmpd1, tmpd2, tsave = p

    # uout = u(t)
    u!(uout, t)

    # hout = d(t-Tau)
    for k=1:length(Tau)
        hout[k] = h(t-Tau[k], k)
    end

    #dx(t) .= A*x(t) + B1*u(t) +B2*d(t-tau)
    mul!(dx, A, x)
    mul!(tmpx, B1, uout)
    dx .+= tmpx
    mul!(tmpx, B2, hout)
    dx .+= tmpx

    return
end

# Save d(t) and y(t) in history
function dde_param_saver(x,t,integrator)
    A, B1, B2, C1, C2, D11, D12, D21, D22, y, Tau, u!, h, uout, hout, tmpx, tmpy1, tmpy2, tmpd1, tmpd2, tsave = integrator.p

    nd = length(Tau)

    # uout = u(t)
    u!(uout, t)

    # hout = d(t-Tau)
    for k=1:length(Tau)
        hout[k] = h(t-Tau[k], k) 
    end

    # tmpd1 = d(t) = C2*x + D21*u(t) + D22*d(t-Tau)
    mul!(tmpd1, C2, x)
    mul!(tmpd2, D21, uout)
    tmpd1 .+= tmpd2
    mul!(tmpd2, D22, hout)
    tmpd1 .+= tmpd2

    # Save delay d(t)
    for k=1:length(Tau)
        h[t,k] = tmpd1[k]
    end

    if t in tsave
        # tmpy1 = y(t) = C1*x + D11*u(t) + D12*d(t-Tau)
        mul!(tmpy1, C1, x)
        mul!(tmpy2, D11, uout)
        tmpy1 .+= tmpy2
        mul!(tmpy2, D12, tmpd2)
        tmpy1 .+= tmpy2
        push!(y, copy(tmpy1))
        println("Saving y t: $t")
    end

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
    tmpx = similar(x0, nx)              # in place storage for x
    tmpy1 = similar(x0, ny)             # in place storage for output
    tmpy2 = similar(x0, ny)             # in place storage for output
    tmpd1 = similar(x0, nd)             # in place storage for delays
    tmpd2 = similar(x0, nd)             # in place storage for delays
    y = Vector{T}[] #Matrix{T}(undef, ny, nt)        # Output matrix
    x = Matrix{T}(undef, nx, nt)        # State matrix

    # History function for delayed variables
    h = HistoryFunction(T, (t,indx) -> 0.0, 0.0, Tau)

    p = (A, B1, B2, C1, C2, D11, D12, D21, D22, y, Tau, u!, h, uout, hout, tmpx, tmpy1, tmpy2, tmpd1, tmpd2, t)
    # This callback computes and stores the delay term
    cb1 = FunctionCallingCallback(dde_param_saver,
                        funcat = t,
                        func_everystep=true,
                        func_start = false)
    # TODO Check what the real limit is on the stepsize
    tau_min = length(Tau)>0 ? minimum(Tau)/2 : 10.0
    tau_min = min(tau_min, dt/100) # TODO This seems to increase the accuracy linearly
    cb2 = StepsizeLimiter((u,p,t) -> tau_min)
    prob = ODEProblem{true}(dde_param, x0,
                (T(t[1]), T(t[end])),
                p,
                callback=CallbackSet(cb1,cb2))#, constant_lags=Tau, stop_at=[0;Tau])

    sol = DelayDiffEq.solve(prob, alg; saveat=t, kwargs..., abstol=1e-10, reltol=1e-10)

    solu = sol.u::Vector{Vector{T}} # the states are labeled u in DelayDiffEq

    for k = 1:nt
        x[:,k] .= sol.u[k]
    end
    # TODO Handle y better
    return reduce(hcat,y[1:end-1])', t, x'
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
