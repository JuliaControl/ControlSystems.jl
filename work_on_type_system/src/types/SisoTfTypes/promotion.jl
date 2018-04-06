Base.promote_rule(::Type{SisoRational{T1}}, ::Type{SisoRational{T2}}) where {T1, T2} = SisoRational{promote_type(T1, T2)}
Base.promote_rule(::Type{SisoZpk{T1}}, ::Type{SisoZpk{T2}}) where {T1, T2} = SisoZpk{promote_type(T1, T2)}

Base.promote_rule(::Type{SisoRational{T1}}, ::Type{SisoZpk{T2}}) where {T1, T2} = SisoRational{promote_type(T1, T2)}

Base.promote_rule(::Type{<:SisoTf}, ::Type{SisoGeneralized}) = SisoGeneralized # ??



Base.promote_rule(::Type{SisoRational{T1}}, ::Type{T2}) where {T1<:Number, T2<:Number} = SisoRational{promote_type(T1, T2)}
Base.promote_rule(::Type{SisoZpk{T1}}, ::Type{T2}) where {T1<:Number, T2<:Number} = SisoZpk{promote_type(T1, T2)}
Base.promote_rule(::Type{SisoGeneralized{T1}}, ::Type{T2}) where {T1<:Number, T2<:Number} = SisoGeneralized{promote_type(T1, T2)}
