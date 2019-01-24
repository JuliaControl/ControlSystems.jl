function freqresp(sys::DelayLtiSystem, ω::AbstractVector{T}) where {T <: Real}
    ny = noutputs(sys)
    nu = ninputs(sys)

    P_fr = ControlSystems.freqresp(sys.P.P, ω);

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



function lsim(sys::DelayLtiSystem, t::AbstractArray{<:Real}; u=(t -> fill(0.0, ninputs(sys))), x0=fill(0.0, nstates(sys)), alg=MethodOfSteps(Tsit5()), kwargs...)
    P = sys.P

    if ~iszero(P.D22)
        error("non-zero D22-matrix block is not supported") # Due to limitations in differential equations
    end

    dt = t[2] - t[1]
    if ~all(diff(t) .≈ dt) # QUESTION Does this work or are there precision problems?
        error("The t-vector should be uniformly spaced, t[2] - t[1] = $dt.") # Perhaps dedicated function for checking this?
    end

    # Slightly more complicated definition since the d signal is not directly available
    dde = function (dx, x, h, p, t)
        dx .= P.A*x + P.B1*u(t)
        for k=1:length(sys.Tau) # Add each of the delayed signals
            dk_delayed = dot(P.C2[k,:], h(p,t-sys.Tau[k])) + dot(P.D21[k,:], u(t-sys.Tau[k]))
            dx .+= P.B2[:, k] * dk_delayed
        end
    end

    h_initial = (p, t) -> zeros(nstates(sys))

    # Saves y (excluding the d(t-τ) contribution) and d
    saved_values = SavedValues(Float64, Tuple{Vector{Float64}, Vector{Float64}})
    cb = SavingCallback((x,t,integrator) -> (P.C1*x + P.D11*u(t), P.C2*x + P.D21*u(t)), saved_values, saveat=t)

    prob = DDEProblem(dde, x0, h_initial, (0.0, 8.0), constant_lags=sys.Tau)

    sol = solve(prob, alg, callback=cb; saveat=t, kwargs...)

    x = sol.u # the states are labeled u in DifferentialEquations
    y = hcat([saved_values.saveval[k][1] for k=1:length(t)]...)
    d = hcat([saved_values.saveval[k][2] for k=1:length(t)]...)

    # Account for the effect of the delayed d-signal on y
    for k=1:length(sys.Tau)
        N_del = Integer(sys.Tau[k] / dt)
        dk = [zeros(N_del); d[k, 1:end-N_del]]

        for j=1:length(t)
            y[:, j] .+= sys.P.D12[:, k] * dk[j]
        end
    end


    t, x, y
end
