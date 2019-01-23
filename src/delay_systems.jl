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

        delay_matrix_fr = Diagonal(exp.(im*sys.Tau*ω[ω_idx])) # Frequency response of the diagonal matrix with delays

        G_fr[ω_idx,:,:] .= P11_fr + P12_fr/(delay_matrix_fr - P22_fr)*P21_fr # The matrix is invertible (?!)
    end

    return G_fr
end



function simulate5(sys::DelayLtiSystem, tspan; u=[], x0=Float64[], alg=MethodOfSteps(Tsit5()), kwargs...)
    P = sys.P

    if ~iszero(P.D22)
        error("non-zero D22-matrix block is not supported") # Due to limitations in differential equations
    end

    nu = ControlSystems.ninputs(sys)
    nx = ControlSystems.nstates(sys)

    u = (u == []) ? t -> fill(0.0, nu) : u
    x0 = (x0 == []) ? fill(0.0, nx) : x0

    dde = function (dx, x, h, p, t)
        dx .= P.A*x + P.B1*u(t)
        for k=1:length(sys.Tau) # Add each of the delayed signals
            dk_delayed = dot(P.C2[k,:], h(p,t-sys.Tau[k])) + dot(P.D21[k,:], u(t-sys.Tau[k]))# + P.D22*(t-Tau[1]
            dx .+= P.B2[:, k] * dk_delayed
        end
        #d_delayed = zeros(size(P.C2,1))
    end

    h_initial = (p, t) -> zeros(nx)


    # Saves y (excluding the d(t-τ) contribution) and d
    saved_values_y = SavedValues(Float64, Tuple{Vector{Float64}, Vector{Float64}})
    cb = SavingCallback((x,t,integrator) -> (P.C1*x + P.D11*u(t), P.C2*x + P.D21*u(t)), saved_values_y, saveat=0:0.02:tspan[2])
#P.C1*x + P.D11*u(t) + P.D12*h(t)
    prob = DDEProblem(dde, x0, h_initial, tspan, constant_lags=sys.Tau)

    sol = solve(prob, alg, callback=cb; saveat=0:0.02:tspan[2], kwargs...)

    t, x = sol.t, sol.u


    t, x, saved_values_y
end
