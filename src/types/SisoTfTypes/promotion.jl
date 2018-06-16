Base.promote_rule(::Type{SisoRational{T1}}, ::Type{SisoRational{T2}}) where {T1, T2} = SisoRational{promote_type(T1, T2)}
Base.promote_rule(::Type{SisoZpk{T1,TR1}}, ::Type{SisoZpk{T2,TR2}}) where {T1, TR1, T2, TR2} = SisoZpk{promote_type(T1, T2), promote_type(TR1,TR2)}

Base.promote_rule(::Type{SisoRational{T1}}, ::Type{SisoZpk{T2}}) where {T1, T2} = SisoRational{promote_type(T1, T2)}


Base.promote_rule(::Type{SisoRational{T1}}, ::Type{T2}) where {T1<:Number, T2<:Number} = SisoRational{promote_type(T1, T2)}
Base.promote_rule(::Type{SisoZpk{T1}}, ::Type{T2}) where {T1<:Number, T2<:Number} = SisoZpk{promote_type(T1, T2)}

Base.promote_op{T<:SisoTf}(::Any, ::Type{T}, ::Type{T}) = T
