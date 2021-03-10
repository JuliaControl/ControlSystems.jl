# """
#     responsetype(T)
# Return the numerical type of the time-response for composite LTISystem type `T`, e.g., `responsetype(ss(1.)) == Float64`
# If `T` is not a type, return `responsetype(typeof(T))`
# """
# responsetype{T<:Real}(::Type{T}) = T
# responsetype{T<:Complex}(::Type{T}) = promote_type(T.types...)
# responsetype(x) = typeof(x) # Fallback method
# responsetype{T}(::Type{T}) = T # Catch all


# Abstract type pyramid =============================================================

# NOTE: a lot that doesnt seem to be used...
#Base.promote_rule(::Type{StateSpace{T1}}, ::Type{StateSpace{T2}}) where {T1,T2} = StateSpace{promote_type(T1, T2)}

# NOTE: Is the below thing correct always?
Base.promote_rule(::Type{StateSpace{TE1,T1}}, ::Type{StateSpace{TE2,T2}}) where {TE1,TE2,T1,T2} =
    StateSpace{promote_type(TE1,TE2), promote_type(T1,T2)}

Base.promote_rule(::Type{StateSpace{TE1,T}}, ::Type{HeteroStateSpace{TE2,AT,BT,CT,DT}}) where {TE1,TE2,T,AT,BT,CT,DT} =
    HeteroStateSpace{promote_type(TE1, TE2),promote_type(Matrix{T},AT),promote_type(Matrix{T},BT),promote_type(Matrix{T},CT),promote_type(Matrix{T},DT)}

Base.promote_rule(::Type{TransferFunction{TE1,S1}}, ::Type{StateSpace{TE2,T2}}) where {TE1,TE2,T1,S1<:SisoTf{T1},T2} =
    StateSpace{promote_type(TE1, TE2), promote_type(T1,T2)}

Base.promote_rule(::Type{TransferFunction{TE1,S1}}, ::Type{DelayLtiSystem{T2}}) where {TE1,T1,S1<:SisoTf{T1},T2} =
    DelayLtiSystem{promote_type(T1,T2)}

Base.promote_rule(::Type{StateSpace{TE1,T1}}, ::Type{DelayLtiSystem{T2}}) where {TE1,T1,T2} =
    DelayLtiSystem{promote_type(T1,T2)}

Base.promote_rule(::Type{TransferFunction{TE1,S1}}, ::Type{TransferFunction{TE2,S2}}) where {TE1, TE2, S1, S2} =
    TransferFunction{promote_type(TE1,TE2),promote_type(S1, S2)}


#function Base.promote_rule{T<:StateSpace,P<:TransferFunction}(::Type{T}, ::Type{P})  #where T <: StateSpace where P <: TransferFunction
#    S = promote_type(primitivereal(P), primitivereal(T)) # TODO: this is not correct for P <: zpk, in which case the StateSpace system gets complex matrices
#    StateSpace{Matrix{S}}
#end


# Promotion of Number ==========================================================

Base.promote_rule(::Type{TransferFunction{TE,S}}, ::Type{T2}) where {TE, S<:SisoTf, T2<:Number} =
    TransferFunction{TE,promote_type(S, T2)}

# TODO Figure out a better way
function Base.promote_rule(::Type{StateSpace{TE, T1}}, ::Type{T2}) where {TE, T1, T2<:Number}
    StateSpace{TE, promote_type(T1, T2)}
end

function Base.promote_rule(::Type{SisoZpk{T1,C2}}, ::Type{T2}) where {T1, C2, T2<:Number}
    GainType = promote_type(T1, T2)
    return SisoZpk{GainType, complex(GainType)}
end
Base.promote_rule(::Type{SisoRational{T1}}, ::Type{T2}) where {T1, T2<:Number} = SisoRational{promote_type(T1,T2)}

# Promotion of Matrix ==========================================================

Base.promote_rule(::Type{TransferFunction{TE, SisoZpk{T1,TR1}}}, ::Type{MT}) where {TE, T1, TR1, T2, MT<:AbstractMatrix{T2}} =
    TransferFunction{TE, SisoZpk{T1, promote_type(TR1, T2)}}

Base.promote_rule(::Type{TransferFunction{TE, SisoRational{T1}}}, ::Type{MT}) where {TE, T1, T2, MT<:AbstractMatrix{T2}} =
    TransferFunction{TE, SisoRational{promote_type(T1, T2)}}

Base.promote_rule(::Type{StateSpace{TE, T1}}, ::Type{MT}) where {TE, T1, MT<:AbstractMatrix} =
    StateSpace{TE, promote_type(T,eltype(MT))}

Base.promote_rule(::Type{DelayLtiSystem{T1,S}}, ::Type{MT}) where {T1, S, MT<:AbstractMatrix} =
    DelayLtiSystem{promote_type(T1, eltype(MT)),S}
