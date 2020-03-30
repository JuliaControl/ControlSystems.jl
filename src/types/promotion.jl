# Base.eltype{T}(::Type{Poly{T}}) = T
# Base.eltype{T}(::Type{TransferFunction{T}}) = T
# Base.eltype{T}(::Type{SisoZpk{T}}) = T
# Base.eltype{T}(::Type{SisoRational{T}}) = T
# Base.eltype{T}(::Type{StateSpace{T}}) = T
#


# """
#     responsetype(T)
# Return the numerical type of the time-response for composite LTISystem type `T`, e.g., `responsetype(ss(1.)) == Float64`
# If `T` is not a type, return `responsetype(typeof(T))`
# """
# responsetype{T<:Real}(::Type{T}) = T
# responsetype{T<:Complex}(::Type{T}) = promote_type(T.types...)
# responsetype(x) = typeof(x) # Fallback method
# responsetype{T}(::Type{T}) = T # Catch all


#matrix_type{T}(::Type{StateSpace{T, M{T}}}) = T

# Abstract type pyramid =============================================================

# NOTE: a lot that doesnt seem to be used...
#Base.promote_rule(::Type{StateSpace{T1}}, ::Type{StateSpace{T2}}) where {T1,T2} = StateSpace{promote_type(T1, T2)}

# NOTE: Is the below thing correct always?
function Base.promote_rule(::Type{StateSpace{T1,MT1}}, ::Type{StateSpace{T2,MT2}}) where {T1,T2,MT1,MT2}
    MT = promote_type(MT1, MT2)
    if isconcretetype(MT) # This sometimes fails and produces a Unionall. This condition can be checked at compile time, checked by f(x,y) = (MT = promote_type(x, y); isconcretetype(MT)); @code_llvm f(Int,Float64)
        StateSpace{promote_type(T1, T2),MT}
    else # If fail, fall back to an ordinary matrix
        FT = promote_type(T1, T2)
        StateSpace{FT,Matrix{FT}}
    end
end

Base.promote_rule(::Type{StateSpace{T,MT}}, ::Type{HeteroStateSpace{AT,BT,CT,DT}}) where {T,MT,AT,BT,CT,DT} = HeteroStateSpace{promote_type(MT,AT),promote_type(MT,BT),promote_type(MT,CT),promote_type(MT,DT)}

function Base.promote_rule(::Type{TransferFunction{S1}}, ::Type{StateSpace{T2,MT}}) where {T1,S1<:SisoTf{T1},T2,MT}
    #promote_type(Base.promote_op(/, T1, T1), T2)
    StateSpace{promote_type(T1,T2), promote_type(Matrix{T1},MT)}
end
# NOTE: Perhaps should try to keep matrix structure?

function Base.promote_rule(::Type{TransferFunction{S1}}, ::Type{DelayLtiSystem{T2}}) where {T1,S1<:SisoTf{T1},T2}
    DelayLtiSystem{promote_type(T1,T2)}
end
function Base.promote_rule(::Type{StateSpace{T1,MT}}, ::Type{DelayLtiSystem{T2}}) where {T1,MT,T2}
    DelayLtiSystem{promote_type(T1,T2)}
end

Base.promote_rule(::Type{TransferFunction{S1}}, ::Type{TransferFunction{S2}}) where {S1, S2} = TransferFunction{promote_type(S1, S2)}
#Base.promote_rule(::Type{SisoTf}, ::Type{TransferFunction}) = TransferFunction
#Base.promote_rule(::Type{SisoZpk}, ::Type{TransferFunction}) = TransferFunction
#Base.promote_rule(::Type{SisoRational}, ::Type{TransferFunction}) = TransferFunction


#function Base.promote_rule{T<:StateSpace,P<:TransferFunction}(::Type{T}, ::Type{P})  #where T <: StateSpace where P <: TransferFunction
#    S = promote_type(primitivereal(P), primitivereal(T)) # TODO: this is not correct for P <: zpk, in which case the StateSpace system gets complex matrices
#    StateSpace{Matrix{S}}
#end


# Promotion of Number ==========================================================

#Base.promote_rule{T1<:Number,T2<:Number}(::Type{TransferFunction{S}}, ::T2) = promote_type(T1, T2)

Base.promote_rule(::Type{TransferFunction{S}}, ::Type{T2}) where {S<:SisoTf, T2<:Number} =
    TransferFunction{promote_type(S, T2)}

# TODO Figure out a better way
function Base.promote_rule(::Type{StateSpace{T1, MT}}, ::Type{T2}) where {T1, MT, T2<:Number}
    NewMT = Base.promote_op(*, MT, T2)
    StateSpace{eltype(NewMT), NewMT}
end

Base.promote_rule(::Type{TransferFunction{SisoZpk{T1,TR1}}}, ::Type{M2}) where {T1, TR1, T2, M2<:AbstractMatrix{T2}} = TransferFunction{SisoZpk{T1, promote_type(TR1, T2)}}

Base.promote_rule(::Type{TransferFunction{SisoRational{T1}}}, ::Type{M2}) where {T1, T2, M2<:AbstractMatrix{T2}} = TransferFunction{SisoRational{promote_type(T1, T2)}}

function Base.promote_rule(::Type{StateSpace{T1, MT1}}, ::Type{MT2}) where {T1, MT1, MT2<:AbstractMatrix}
    MT = promote_type(MT1, MT2)
    StateSpace{eltype(MT), MT}
end

Base.promote_rule(::Type{DelayLtiSystem{T1,S}}, ::Type{MT1}) where {T1, S, MT1<:AbstractMatrix} =
    DelayLtiSystem{promote_type(T1, eltype(MT1)),S}

#Base.promote_rule{S<:TransferFunction{<:SisoTf}}(::Type{S}, ::Type{<:Real}) = S

# We want this, but not possible, so hardcode for SisoTypes
#Base.promote_rule(::Type{S}, ::Type{T2}) where {T1, S<:SisoTf{T1}, T2<:Number} = S{promote_type(T1, T2)}

function Base.promote_rule(::Type{SisoZpk{T1,C2}}, ::Type{T2}) where {T1, C2, T2<:Number}
    GainType = promote_type(T1, T2)
    return SisoZpk{GainType, complex(GainType)}
end
Base.promote_rule(::Type{SisoRational{T1}}, ::Type{T2}) where {T1, T2<:Number} = SisoRational{promote_type(T1,T2)}
#Base.promote_rule(::Type{StateSpace{T1, MT}}, ::Type{T2}) where {T1, T2<:Number} = SisoRational{promote_type(T1,T2)}

#Base.promote_rule{T<:SisoZpk}(::Type{T}, ::Type{<:Real}) = T
#Base.promote_rule{T}(::Type{SisoRational{T}}, ::Type{<:Real}) = SisoRational{T}
#Base.promote_rule{T}(::Type{StateSpace{T}}, ::Type{<:Real}) = StateSpace{T}
#Base.promote_rule{S<:TransferFunction{<:SisoTf},T<:Real}(::Type{S}, ::Union{Type{Array{T,2}},Type{Array{T,1}}}) = S


# Less abstract promotions
#Base.promote_rule(::Type{TransferFunction{SisoRational}}, ::Type{TransferFunction{SisoZpk}}) = TransferFunction{SisoZpk} # NOTE: Is this what we want?




# function Base.promote_rule{T1c<:SisoZpk,T2c<:SisoRational}(::Type{T1c}, ::Type{T2c})
#     upper_type = SisoZpk
#     inner_type = promote_type(arraytype(T1c), arraytype(T2c))
#     upper_type{inner_type}
# end



# @show siso_type = promote_type(eltype.(systems)...)
# if isleaftype(siso_type)
#     siso_type = siso_type.name.wrapper
# end
# # array_type = promote_type(arraytype.(systems)...)
# # mat = hcat([(e->convert(siso_type{array_type}, e)).(s.matrix) for s in systems]...)
# @show promote_type([s.matrix for s in systems]...)




#

#
#
# # Promote_op types # NOTE: Needed?
# Base.promote_op{T<:SisoTf}(::Any, ::Type{T}, ::Type{T}) = T
#
# #Just default SisoTf to SisoRational
# SisoTf(args...) = SisoRational(args...)
#
# Base.zero(::Type{<:SisoTf}) = zero(SisoRational)
# Base.zero(::SisoTf) = zero(SisoRational)
# Base.zero(::Type{<:SisoZpk}) = SisoZpk(Float64[],Float64[],0.0)
# Base.zero(::SisoZpk) = Base.zero(SisoZpk)
# Base.zero{T}(::Type{SisoRational{T}}) = SisoRational(zero(Poly{T}), one(Poly{T})) # FIXME: Is the in analogy with how zero for SisoRational is createdzer?
# Base.zero{T}(::SisoRational{T}) = Base.zero(SisoRational{T})
