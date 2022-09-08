function Base.convert(::Type{SisoZpk{T,TR}}, f::SisoRational{T2}) where {T<:Number, TR<:Number, T2<:Number}
    if length(f.num) == 0
        return SisoZpk(T[],TR[],0)
    elseif all(f.den == zero(f.den))
        error("Zero denominator, this should not be possible") # NOTE: Should we keep this?
    end

    return SisoZpk{T,TR}(roots(f.num), roots(f.den), f.num[end]/f.den[end]) # NOTE: polynomials are stored with increasing coefficient orders
end

function Base.convert(::Type{<:SisoZpk}, f::SisoRational{T}) where {T<:Number}
    T1 = Base.promote_op(/, T, T)
    TR = complex(T1) # Type of roots(f.z)
    Base.convert(SisoZpk{T1,TR}, f)
end


function Base.convert(::Type{SisoZpk{T,TR}}, f::SisoZpk) where {T<:Number, TR<:Number}
    return SisoZpk{T,TR}(f.z, f.p, f.k)
end


function Base.convert(::Type{SisoRational{T}}, f::SisoRational) where {T<:Number}
    return SisoRational{T}(convert(Polynomial{T}, f.num), convert(Polynomial{T}, f.den))
end

function Base.convert(::Type{SisoRational{T1}}, f::SisoZpk{T2,TR}) where {T1<:Number, T2<:Number,TR<:Number}
    if T2 <: Real
        num = prod(roots2real_poly_factors(f.z))*f.k
        den = prod(roots2real_poly_factors(f.p))
    else
        num = prod(roots2poly_factors(f.z))*f.k
        den = prod(roots2poly_factors(f.p))
    end
    return SisoRational{T1}(num, den)
end

function Base.convert(::Type{SisoRational}, f::SisoZpk{T1,TR}) where {T1<:Number, TR<:Number}
    if T1 <: Real
        T = promote_type(T1,typeof(real(one(TR))))
        # Alternative version if T1 rich enough
        # T = T1
        convert(SisoRational{T}, f)
    else
        convert(SisoRational{promote_type(T1,TR)}, f)
    end
end


Base.convert(::Type{S}, b::T2) where {T1, S <: SisoTf{T1}, T2<:Number} = convert(T1, b)*one(S)
