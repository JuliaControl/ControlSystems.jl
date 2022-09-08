Base.promote_rule(::Type{SisoRational{T1}}, ::Type{SisoRational{T2}}) where {T1, T2} = SisoRational{promote_type(T1, T2)}
Base.promote_rule(::Type{SisoZpk{T1,TR1}}, ::Type{SisoZpk{T2,TR2}}) where {T1, TR1, T2, TR2} = SisoZpk{promote_type(T1, T2), promote_type(TR1,TR2)}

function Base.promote_rule(::Type{SisoRational{T1}}, ::Type{SisoZpk{T2,TR2}}) where {T1, T2, TR2<:Number}
    Tnew = promote_type(T1, T2)
    TRnew = promote_type(TR2, complex(Tnew))
    SisoZpk{Tnew, TRnew}
end


Base.promote_rule(::Type{SisoRational{T1}}, ::Type{T2}) where {T1<:Number, T2<:Number} = SisoRational{promote_type(T1, T2)}

function Base.promote_rule(::Type{SisoZpk{T1,TR1}}, ::Type{T2}) where {T1<:Number, TR1<:Number, T2<:Number}
    Tnew = promote_type(T1, T2)
    TRnew = promote_type(TR1, complex(Tnew))
    # TODO Not obvious that we want to promote to complex poles?
    SisoZpk{Tnew, TRnew}
end
