function freqresp(sys::DelayLtiSystem, ω::AbstractVector{T}) where {T <: Real}
    ny = noutputs(sys)
    nu = ninputs(sys)

    P_fr = ControlSystems.freqresp(sys.P, ω);

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
