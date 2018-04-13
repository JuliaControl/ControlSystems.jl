function Base.convert(::Type{SisoZpk{T,TR}}, f::SisoRational) where {T<:Number, TR<:Number}
    if length(f.num) == 0
        return SisoZpk(TR[],TR[],0)
    elseif all(f.den == zero(f.den))
        error("Zero denominator, this should not be possible") # NOTE: Should we keep this?
    end

    return SisoZpk{T,TR}(roots(f.num), roots(f.den), f.num[end]/f.den[end]) # NOTE: polynomials are stored with increasing coefficient orders
end
function Base.convert(::Type{<:SisoZpk}, f::SisoRational)
    T = numeric_type(f)
    TR = promote_type(T, Complex128) # Borde vara typen för roots(f.z)
    Base.convert(SisoZpk{T,TR}, f)
end


function Base.convert(::Type{SisoZpk{T,TR}}, f::SisoZpk) where {T<:Number, TR<:Number}
    return SisoZpk{T,TR}(f.z, f.p, f.k)
end


function Base.convert(::Type{<:SisoRational{T}}, f::SisoRational) where T
    return SisoRational{T}(convert(Poly{T}, f.num), convert(Poly{T}, f.den))
end

function Base.convert(::Type{<:SisoRational}, f::SisoZpk)
    T = numeric_type(f)
    if isreal(one(T)) # FIXME: This probably doesn't do the trick
        num = prod(roots2real_poly_factors(f.z))*f.k
        den = prod(roots2real_poly_factors(f.p))
    else
        num = prod(roots2complex_poly_factors(f.z))*f.k
        den = prod(roots2complex_poly_factors(f.p))
    end

    return SisoRational{T}(num, den)
end



Base.convert(::Type{S}, b) where {T, S <: SisoTf{T}} = convert(b, T)*one(S)
