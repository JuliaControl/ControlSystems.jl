## User should just use TransferFunction
struct SisoRational{VT<:AbstractVector{<: Number}} <: SisoTf
    num::Poly{VT}
    den::Poly{VT}
    function SisoRational(num::Poly{VT}, den::Poly{VT}) where VT <: AbstractVector{<: Number}
        if all(den == zero(den))
            error("Zero denominator")
        elseif all(num == zero(num))
            # The numerator is zero, make the denominator 1
            den = one(den)
        end
        new{VT}(num, den)
    end
end
SisoRational(num::Vector{T}, den::Vector{T}) where T <: Number = SisoRational(Poly(num), Poly(den))
function SisoRational(num::Poly{T1}, den::Poly{T2}) where T1 <: AbstractVector{<: Number} where T2 <: AbstractVector{<: Number}
    T = promote_type(eltype(T1),eltype(T2))
    SisoRational(Poly(map(T,num.a)), Poly(map(T,den.a)))
end

function minreal(sys::SisoRational, eps::Real=sqrt(eps()))
    return SisoRational(minreal(SisoZpk(sys), eps))
end

function print_siso(io::IO, t::SisoRational, var=:s)
    # Convert the numerator and denominator to strings
    numstr = sprint(print_poly, t.num, var)
    denstr = sprint(print_poly, t.den, var)

    # Figure out the length of the separating line
    len_num = length(numstr)
    len_den = length(denstr)
    dashcount = max(len_num, len_den)

    # Center the numerator or denominator
    if len_num < dashcount
        numstr = "$(repeat(" ", div(dashcount - len_num, 2)))$numstr"
    else
        denstr = "$(repeat(" ", div(dashcount - len_den, 2)))$denstr"
    end
    println(io, numstr)
    println(io, repeat("-", dashcount))
    println(io, denstr)
end

function print_compact(io::Base.IO, t::SisoRational, var=:s)
    numstr = sprint(print_poly, t.num, var)
    denstr = sprint(print_poly, t.den, var)
    println(io, "($numstr)/($denstr)")
end

Base.promote_rule{T<:Real}(::Type{SisoRational}, ::Type{T}) = SisoRational
Base.convert(::Type{SisoRational}, b::Real) = SisoRational([b], [1])

Base.zero(::Type{SisoRational}) = SisoRational(zero(Poly{Vector{Float64}}), one(Poly{Vector{Float64}}))
Base.zero(::SisoRational) = Base.zero(SisoRational)

Base.length(t::SisoRational) = max(length(t.num), length(t.den))

function Base.num(t::SisoRational)
    lt = length(t)
    n = zeros(eltype(t.num), lt)
    n[(lt - length(t.num) + 1):end] = t.num[:]
    return n
end

function Base.den(t::SisoRational)
    lt = length(t)
    d = zeros(eltype(t.den), lt)
    d[(lt - length(t.den) + 1):end] = t.den[:]
    return d
end

denpoly(G::SisoRational) = Poly(den(G))
numpoly(G::SisoRational) = Poly(num(G))

function evalfr(sys::SisoRational, s::Number)
    S = promote_type(typeof(s), Float64)
    den = polyval(sys.den, s)
    if den == zero(S)
        convert(S, Inf)
    else
        polyval(sys.num, s)/den
    end
end

tzero(t::SisoRational) = roots(t.num)
pole(t::SisoRational) = roots(t.den)

==(t1::SisoRational, t2::SisoRational) = (t1.num == t2.num && t1.den == t2.den)
# We might want to consider alowing scaled num and den as equal
function isapprox(t1::SisoRational, t2::SisoRational; rtol::Real=sqrt(eps()), atol::Real=0)
    isapprox(t1.num,t2.num, rtol=rtol, atol=atol) && isapprox(t1.den, t2.den, rtol=rtol, atol=atol)
end

+(t1::SisoRational, t2::SisoRational) = SisoRational(t1.num*t2.den + t2.num*t1.den, t1.den*t2.den)
+(t::SisoRational, n::Real) = SisoRational(t.num + n*t.den, t.den)
+(n::Real, t::SisoRational) = t + n
#.+(t::SisoRational, n::Real) = t + n
#.+(n::Real, t::SisoRational) = t + n

-(t1::SisoRational, t2::SisoRational) = SisoRational(t1.num*t2.den - t2.num*t1.den, t1.den*t2.den)
-(n::Real, t::SisoRational) = SisoRational(n*t.den - t.num, t.den)
-(t::SisoRational, n::Real) = +(t, -n)
#.-(t::SisoRational, n::Real) = -(t, n)
#.-(n::Real, t::SisoRational) = -(n, t)

-(t::SisoRational) = SisoRational(-t.num, t.den)

*(t1::SisoRational, t2::SisoRational) = SisoRational(t1.num*t2.num, t1.den*t2.den)
*(t::SisoRational, n::Real) = SisoRational(t.num*n, t.den)
*(n::Real, t::SisoRational) = *(t, n)
#.*(t1::SisoRational, t2::SisoRational) = *(t1, t2)
#.*(t::SisoRational, n::Real) = *(t, n)
#.*(n::Real, t::SisoRational) = *(t, n)

/(n::Real, t::SisoRational) = SisoRational(n*t.den, t.num)
/(t::SisoRational, n::Real) = t*(1/n)
/(t1::SisoRational, t2::SisoRational) = t1*(1/t2)
