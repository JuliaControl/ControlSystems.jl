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

        delay_matrix_inv_fr = Diagonal(exp.(im*sys.Tau*ω[ω_idx])) # Frequency response of the diagonal matrix with delays
        # Inverse of the delay matrix, so there should not be any minus signs in the exponents

        G_fr[ω_idx,:,:] .= P11_fr + P12_fr/(delay_matrix_inv_fr - P22_fr)*P21_fr # The matrix is invertible (?!)
    end

    return G_fr
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
function lsim(sys::DelayLtiSystem{T,S}, u, t::AbstractArray{<:Real}; x0=fill(zero(T), nstates(sys)), alg=MethodOfSteps(Tsit5()), kwargs...) where {T,S}
    # Make u! in-place function of u
    u! = if isa(u, Number) || isa(u,AbstractVector) # Allow for u to be a constant number or vector
        (uout, t) -> uout .= u
    elseif DiffEqBase.isinplace(u, 2)               # If u is an inplace (more than 1 argument function)
        u
    else                                            # If u is a regular u(t) function
        (out, t) -> (out .= u(t))
    end
    _lsim(sys, u!, t, x0, alg; kwargs...)
end

function dde_param(dx, x, h!, p, t)
    A, B1, B2, C2, D21, Tau, u!, uout, hout, tmp = p

    u!(uout, t)     # uout = u(t)

    #dx .= A*x + B1*ut
    mul!(dx, A, x)
    mul!(tmp, B1, uout)
    dx .+= tmp

    @views for k=1:length(Tau)     # Add each of the delayed signals
        u!(uout, t-Tau[k])      # uout = u(t-tau[k])
        h!(hout, p, t-Tau[k])
        dk_delayed = dot(C2[k,:], hout) + dot(D21[k,:], uout)
        dx .+= B2[:, k] .* dk_delayed
    end
    return
end

# TODO Discontinuities in u are not handled well yet.
function _lsim(sys::DelayLtiSystem{T,S}, Base.@nospecialize(u!), t::AbstractArray{<:Real}, x0::Vector{T}, alg; kwargs...) where {T,S}
    P = sys.P

    if ~iszero(P.D22)
        throw(ArgumentError("non-zero D22-matrix block is not supported")) # Due to limitations in differential equations
    end

    t0 = first(t)
    dt = t[2] - t[1]
    if ~all(diff(t) .≈ dt) # QUESTION Does this work or are there precision problems?
        error("The t-vector should be uniformly spaced, t[2] - t[1] = $dt.") # Perhaps dedicated function for checking this?
    end

    # Get all matrices to save on allocations
    A, B1, B2, C1, C2, D11, D12, D21, D22 = P.A, P.B1, P.B2, P.C1, P.C2, P.D11, P.D12, P.D21, P.D22
    Tau = sys.Tau

    hout = fill(zero(T), nstates(sys))  # in place storage for h
    uout = fill(zero(T), ninputs(sys)) # in place storage for u
    tmp = similar(x0)

    h!(out, p, t) = (out .= 0)      # History function

    p = (A, B1, B2, C2, D21, Tau, u!, uout, hout, tmp)
    prob = DDEProblem{true}(dde_param, x0, h!, (T(t[1]), T(t[end])), p, constant_lags=sys.Tau)

    sol = DelayDiffEq.solve(prob, alg; saveat=t, kwargs...)

    x = sol.u::Vector{Vector{T}} # the states are labeled u in DelayDiffEq

    y = Matrix{T}(undef, noutputs(sys), length(t))
    d = Matrix{T}(undef, size(C2,1), length(t))
    # Build y signal (without d term)
    for k = 1:length(t)
        u!(uout, t[k])
        y[:,k] = C1*x[k] + D11*uout
        #d[:,k] = C2*x[k] + D21*uout
    end

    dtmp = Vector{T}(undef, size(C2,1))
    # Function to evaluate d(t)_i at an arbitrary time
    # X is continuous, so interpoate, u is not
    function dfunc!(tmp::Array{T1}, t, i) where T1
        tmp .= if t < t0
            T1(0)
        else
            sol(t) # solution object has built in interpolator of same order as solver.
        end
        u!(uout, t)
        return dot(view(C2,i,:),tmp) + dot(view(D21,i,:),uout)
    end

    # Account for the effect of the delayed d-signal on y
    for k=1:length(Tau)
        for j=1:length(t)
            di = dfunc!(tmp, t[j] - Tau[k], k)
            y[:, j] .+= view(D12,:, k) .* di
        end
    end

    return y', t, reduce(hcat ,x)'
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
function _default_Ts(sys::DelayLtiSystem)
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

iscontinuous(sys::DelayLtiSystem) = true

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
