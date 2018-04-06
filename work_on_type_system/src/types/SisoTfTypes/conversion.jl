function Base.convert(::Type{SisoZpk{T,TR}}, f::SisoRational) where {T<:Number, TR<:Number}
    if length(f.num) == 0
        return SisoZpk(TR[],TR[],0)
    elseif all(f.den == zero(f.den))
        error("Zero denominator, this should not be possible") # NOTE: Should we keep this?
    end

    return SisoZpk{T,TR}(roots(f.num), roots(f.den), f.num[1]/f.den[1])
end
function Base.convert(::Type{<:SisoZpk}, f::SisoRational)
    T = numbertype(typeof(f))
    TR = promote_type(T, Complex128) # Borde vara typen för roots(f.z)
    Base.convert(SisoZpk{T,TR}, f)
end



function Base.convert(::Type{<:SisoRational{T}}, f::SisoRational) where T
    return SisoRational{T}(convert(Poly{T}, f.num), convert(Poly{T}, f.den))
end


function Base.convert(::Type{<:SisoRational}, f::SisoZpk)
    T = numbertype(typeof(f))
    if isreal(T)
        num = prod(roots2real_poly_factors(f.z))*f.k
        den = prod(roots2real_poly_factors(sys.p))
    else
        error("Not implemented yet")
    end

    return SisoRational{T}(num, den)
end



function Base.convert(::Type{<:SisoRational}, sys::SisoZpk)
    num = prod(zp2polys(sys.z))*sys.k
    den = prod(zp2polys(sys.p))
    return SisoRational(num, den)
end



Base.convert(::Type{S}, b) where {T, S <: SisoTf{T}} = convert(b, T)*one(S)
