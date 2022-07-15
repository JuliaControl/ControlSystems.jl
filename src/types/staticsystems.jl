const StaticStateSpace{TE, ny, nu, nx} = HeteroStateSpace{TE,
    <:SMatrix{nx,nx},
    <:SMatrix{nx,nu},
    <:SMatrix{ny,nx},
    <:SMatrix{ny,nu}
}


to_static(a::Number) = a
to_static(a::AbstractArray) = SMatrix{size(a)...}(a)
to_static(a::SArray) = a
to_sized(a::Number) = a
to_sized(a::AbstractArray) = SizedArray{Tuple{size(a)...}}(a)

function HeteroStateSpace(A,B,C,D,Ts=0,f::F=to_static) where F
    HeteroStateSpace(f(A),f(B),f(C),f(D),Ts)
end

HeteroStateSpace(s,f) = HeteroStateSpace(s.A,s.B,s.C,s.D,s.timeevol,f)

ControlSystems._string_mat_with_headers(a::SizedArray) = ControlSystems._string_mat_with_headers(Matrix(a))



"""
    StaticStateSpace(sys::AbstractStateSpace)

Return a [`HeteroStateSpace`](@ref) where the system matrices are of type SMatrix.

*NOTE: This function is fundamentally type unstable.* For maximum performance, create the static system manually, or make use of the function-barrier technique.
"""
StaticStateSpace(sys::AbstractStateSpace) = HeteroStateSpace(sys, to_static)

"""
    to_sized(sys::AbstractStateSpace)

Return a [`HeteroStateSpace`](@ref) where the system matrices are of type SizedMatrix.

*NOTE: This function is fundamentally type unstable.* For maximum performance, create the sized system manually, or make use of the function-barrier technique.
"""
to_sized(sys::AbstractStateSpace) = HeteroStateSpace(sys, to_sized)

StaticStateSpace(G::TransferFunction) = StaticStateSpace(ss(G))

# function to_static(sys::DelayLtiSystem)
#     innerP = to_static(sys.P.P)
#     partP = PartitionedStateSpace(innerP, sys.P.nu1, sys.P.ny1)
#     DelayLtiSystem(partP, sys.Tau)
# end


## Feedback
# Special method to make sure we handle static systems in a type stable way


function HeteroStateSpace(D::SArray, Ts=nothing)
    HeteroStateSpace(D, Ts === nothing ? Continuous() : Discrete(Ts))
end

function HeteroStateSpace(D::SArray, timeevol::TimeEvolution)
    HeteroStateSpace(
        SMatrix{0,0,eltype(D),0}(),
        SMatrix{0,size(D,2),eltype(D),0}(),
        SMatrix{size(D,1),0,eltype(D),0}(),
        D,
        timeevol
    )
end

function Base.promote_rule(::Type{N}, ::Type{<:HeteroStateSpace{TE}}) where {N <: Number, TE}
    HeteroStateSpace{TE}
end

function Base.convert(::Type{<:HeteroStateSpace{TE}}, n::Number) where TE
    if TE <: Continuous
        HeteroStateSpace(SMatrix{1,1,typeof(n),1}(n), Continuous())
    else
        HeteroStateSpace(SMatrix{1,1,typeof(n),1}(n), Discrete(UNDEF_SAMPLEPETIME))
    end
end

function Base.promote_rule(::Type{StateSpace{TE1,T}}, ::Type{HeteroStateSpace{TE2,AT,BT,CT,DT}}) where {TE1,TE2,T,AT<:SMatrix,BT<:SMatrix,CT<:SMatrix,DT<:SMatrix}
    StateSpace{promote_type(TE1, TE2), promote_type(T, eltype(AT), eltype(BT), eltype(CT), eltype(DT))}
end


function Base.promote_rule(::Type{TransferFunction{TE1,T}}, ::Type{HeteroStateSpace{TE2,AT,BT,CT,DT}}) where {TE1,TE2,T,AT<:SMatrix,BT<:SMatrix,CT<:SMatrix,DT<:SMatrix}
    StateSpace{promote_type(TE1, TE2), promote_type(numeric_type(T), eltype(AT), eltype(BT), eltype(CT), eltype(DT))}
end


function feedback(sys1::StaticStateSpace, sys2::StaticStateSpace;
    U1=:, Y1=:, U2=:, Y2=:, W1=:, Z1=:, W2=SVector{0,Int}(), Z2=SVector{0,Int}(),
    Wperm=:, Zperm=:, pos_feedback::Bool=false)

    timeevol = common_timeevol(sys1,sys2)

    if !(isa(Y1, Colon) || allunique(Y1)); @warn "Connecting single output to multiple inputs Y1=$Y1"; end
    if !(isa(Y2, Colon) || allunique(Y2)); @warn "Connecting single output to multiple inputs Y2=$Y2"; end
    if !(isa(U1, Colon) || allunique(U1)); @warn "Connecting multiple outputs to a single input U1=$U1"; end
    if !(isa(U2, Colon) || allunique(U2)); @warn "Connecting a single output to multiple inputs U2=$U2"; end

    if (U1 isa Colon ? size(sys1, 2) : length(U1)) != (Y2 isa Colon ? size(sys2, 1) : length(Y2))
        error("Lengths of U1 ($U1) and Y2 ($Y2) must be equal")
    end
    if (U2 isa Colon ? size(sys2, 2) : length(U2)) != (Y1 isa Colon ? size(sys1, 1) : length(Y1))
        error("Lengths of U1 ($U2) and Y2 ($Y1) must be equal")
    end

    α = pos_feedback ? 1 : -1 # The sign of feedback

    s1_B1 = sys1.B[:,W1]
    s1_B2 = sys1.B[:,U1]
    s1_C1 = sys1.C[Z1,:]
    s1_C2 = sys1.C[Y1,:]
    s1_D11 = sys1.D[Z1,W1]
    s1_D12 = sys1.D[Z1,U1]
    s1_D21 = sys1.D[Y1,W1]
    s1_D22 = sys1.D[Y1,U1]

    s2_B1 = sys2.B[:,W2]
    s2_B2 = sys2.B[:,U2]
    s2_C1 = sys2.C[Z2,:]
    s2_C2 = sys2.C[Y2,:]
    s2_D11 = sys2.D[Z2,W2]
    s2_D12 = sys2.D[Z2,U2]
    s2_D21 = sys2.D[Y2,W2]
    s2_D22 = sys2.D[Y2,U2]

    if iszero(s1_D22) || iszero(s2_D22)
        A = [[sys1.A + α*s1_B2*s2_D22*s1_C2        α*s1_B2*s2_C2];
                 [s2_B2*s1_C2            sys2.A + α*s2_B2*s1_D22*s2_C2]]

        B = [[s1_B1 + α*s1_B2*s2_D22*s1_D21        α*s1_B2*s2_D21];
                      [s2_B2*s1_D21            s2_B1 + α*s2_B2*s1_D22*s2_D21]]
        C = [[s1_C1 + α*s1_D12*s2_D22*s1_C2        α*s1_D12*s2_C2];
                      [s2_D12*s1_C2           s2_C1 + α*s2_D12*s1_D22*s2_C2]]
        D = [[s1_D11 + α*s1_D12*s2_D22*s1_D21        α*s1_D12*s2_D21];
                      [s2_D12*s1_D21           s2_D11 + α*s2_D12*s1_D22*s2_D21]]
    else
        R1 = try
            inv(α*I - s2_D22*s1_D22) # slightly faster than α*inv(I - α*s2_D22*s1_D22)
        catch
            error("Ill-posed feedback interconnection,  I - α*s2_D22*s1_D22 or I - α*s2_D22*s1_D22 not invertible")
        end

        R2 = try
            inv(I - α*s1_D22*s2_D22)
        catch
            error("Ill-posed feedback interconnection,  I - α*s2_D22*s1_D22 or I - α*s2_D22*s1_D22 not invertible")
        end

        A = [[sys1.A + s1_B2*R1*s2_D22*s1_C2        s1_B2*R1*s2_C2];
                 [s2_B2*R2*s1_C2            sys2.A + α*s2_B2*R2*s1_D22*s2_C2]]

        B = [[s1_B1 + s1_B2*R1*s2_D22*s1_D21        s1_B2*R1*s2_D21];
                     [s2_B2*R2*s1_D21            s2_B1 + α*s2_B2*R2*s1_D22*s2_D21]]
        C = [[s1_C1 + s1_D12*R1*s2_D22*s1_C2        s1_D12*R1*s2_C2];
                     [s2_D12*R2*s1_C2           s2_C1 + α*s2_D12*R2*s1_D22*s2_C2]]
        D = [[s1_D11 + s1_D12*R1*s2_D22*s1_D21        s1_D12*R1*s2_D21];
                     [s2_D12*R2*s1_D21           s2_D11 + α*s2_D12*R2*s1_D22*s2_D21]]
    end

    Dfinal = D[Zperm, Wperm]
    return HeteroStateSpace(A, B[:, Wperm], C[Zperm,:], Dfinal, timeevol)

end

function *(sys1::StaticStateSpace, sys2::StaticStateSpace)
    #Check dimension alignment
    #Note: sys1*sys2 = y <- sys1 <- sys2 <- u
    if (sys1.nu != sys2.ny) && (sys1.nu == 1 || sys2.ny == 1)
        throw(DimensionMismatch("sys1*sys2: sys1 must have same number of inputs as sys2 has outputs. If you want to broadcast a scalar system to a diagonal system, use broadcasted multiplication sys1 .* sys2"))
    end
    sys1.nu == sys2.ny || throw(DimensionMismatch("sys1*sys2: sys1 must have same number of inputs as sys2 has outputs"))
    timeevol = common_timeevol(sys1,sys2)
    T = promote_type(numeric_type(sys1), numeric_type(sys2))

    O = @SMatrix zeros(size(sys2.A, 1), size(sys1.A, 2))
    A = [[sys1.A    sys1.B*sys2.C];
         [O  sys2.A]]
    B = [sys1.B*sys2.D ; sys2.B]
    C = [sys1.C   sys1.D*sys2.C]
    D = sys1.D*sys2.D
    return HeteroStateSpace(A, B, C, D, timeevol)
end


function *(sys1::StaticStateSpace, sys2::AbstractStateSpace)
    sys1*StaticStateSpace(sys2)
end

function *(sys1::AbstractStateSpace, sys2::StaticStateSpace)
    StaticStateSpace(sys1)*sys2
end



@autovec () function freqresp_nohess!(R::Array{T,3}, sys::StaticStateSpace, w_vec::AbstractVector{W}) where {T, W <: Real}
    ny, nu = size(sys)
    @boundscheck size(R) == (ny,nu,length(w_vec))
    nx = sys.nx
    if nx == 0 # Only D-matrix
        @inbounds for i in eachindex(w_vec)
            R[:,:,i] .= sys.D
        end
        return R
    end
    A,B,C,D = ssdata(sys)
    te = sys.timeevol
    
    let R=R, A=A, B=B, C=C, D=D, te=te
        @inbounds Polyester.@batch for i in eachindex(w_vec)
            Ri = @views R[:,:,i]
            copyto!(Ri,D) # start with the D-matrix
            isinf(w_vec[i]) && continue
            w = _freq(w_vec[i], te)
            Ac = A - w*I
            Bc = Ac \ B # Bc = (A - w*I)\B 
            Ri .-= C*Bc # - rather than + since (A - w*I) instead of (w*I - A)
        end
    end
    R
end
