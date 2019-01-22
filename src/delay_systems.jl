function freqresp(sys::DelayLtiSystem, ω::AbstractVector{T}) where {T <: Real}
    ny = noutputs(sys)
    nu = ninputs(sys)

    P_fr = ControlSystems.freqresp(sys.P, ω);

    println(ω)

    # FIXME: Different dimensions compared to standard freqresp
    G_fr = zeros(eltype(P_fr), ny, nu, length(ω))

    for ω_idx=1:length(ω)
        P11_fr = P_fr[ω_idx, 1:ny, 1:nu]
        P12_fr = P_fr[ω_idx, 1:ny, nu+1:end]
        P21_fr = P_fr[ω_idx, ny+1:end, 1:nu]
        P22_fr = P_fr[ω_idx, ny+1:end, nu+1:end]

        # FIXME: when there is no delays...

        delay_vect_fr = Base.exp.(im*sys.Tau*ω[ω_idx]) # Frequency response of the block diagonal matrix

        G_fr[:,:,ω_idx] .= P11_fr + P12_fr/(Diagonal(delay_vect_fr) - P22_fr)*P21_fr # The matrix is invertible (?!)
    end

    return G_fr
end
