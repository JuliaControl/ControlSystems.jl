## User should just use TransferFunction
struct SisoRational{T} <: SisoTf{T}
    num::Poly{T}
    den::Poly{T}
    function SisoRational{T}(num::Poly{T}, den::Poly{T}) where T <: Number
        if all(den == zero(den))
            error("Cannot create SisoRational with zero denominator")
        elseif all(num == zero(num))
            # The numerator is zero, make the denominator 1
            den = one(den)
        end
        new{T}(num, den)
    end
end
function SisoRational(num::Poly{T1}, den::Poly{T2}) where T1 <: Number where T2 <: Number
    T = promote_type(T1,T2)
    SisoRational{T}(Poly{T}(num.a), Poly{T}(den.a))
end
function SisoRational{T}(num::AbstractVector, den::AbstractVector) where T <: Number # NOTE: Typearguemnts on the parameters?
    SisoRational{T}(Poly{T}(reverse(num)), Poly{T}(reverse(den)))
end
function SisoRational(num::AbstractVector{T1}, den::AbstractVector{T2}) where T1 <: Number where T2 <: Number
    T = promote_type(T1,T2)
    SisoRational{T}(num, den)
end
# NOTE: How many of these above are actually needed?
# TODO: Add method for scalar inputs


Base.zero(::Type{SisoRational{T}}) where T = SisoRational{T}([zero(T)], [one(T)])
Base.one(::Type{SisoRational{T}}) where T = SisoRational{T}([one(T)], [one(T)])

Base.one(f::SisoRational) = one(typeof(f))
Base.zero(f::SisoRational) = zero(typeof(f))

isproper(f::SisoRational) = (length(f.num) <= length(f.den))

function minreal(sys::SisoRational, eps::Real=sqrt(eps()))
    return SisoRational(minreal(SisoZpk(sys), eps))
end

function print_siso(io::IO, f::SisoRational, var=:s)
    # Convert the numerator and denominator to strings
    numstr = sprint(printpolyfun(var), f.num)
    denstr = sprint(printpolyfun(var), f.den)

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

# TODO is this working?
function print_compact(io::Base.IO, f::SisoRational, var)
    numstr = sprint(print_poly, f.num)
    denstr = sprint(print_poly, f.den)
    println(io, "($numstr)/($denstr)")
end


Base.num(f::SisoRational) = reverse(coeffs(f.num))
Base.den(f::SisoRational) = reverse(coeffs(f.den))

denpoly(f::SisoRational) = f.den
numpoly(f::SisoRational) = f.num

tzero(f::SisoRational) = roots(f.num)
pole(f::SisoRational) = roots(f.den)

function evalfr(f::SisoRational, s::Number)
    S = promote_type(typeof(s), Float64)
    den = polyval(f.den, s)
    if den == zero(S)
        convert(S, Inf)
    else
        polyval(f.num, s)/den
    end
end

==(f1::SisoRational, f2::SisoRational) = (f1.num == f2.num && f1.den == f2.den) # NOTE: Not in analogy with how it's done for SisoZpk

# We might want to consider alowing scaled num and den as equal
function isapprox(f1::SisoRational, f2::SisoRational; rtol::Real=sqrt(eps()), atol::Real=0)
    isapprox(f1.num,f2.num, rtol=rtol, atol=atol) && isapprox(f1.den, f2.den, rtol=rtol, atol=atol)
end

+(f1::SisoRational, f2::SisoRational) = SisoRational(f1.num*f2.den + f2.num*f1.den, f1.den*f2.den)
+(f::SisoRational, n::Real) = SisoRational(f.num + n*f.den, f.den)
+(n::Real, f::SisoRational) = f + n
#.+(f::SisoRational, n::Real) = t + n
#.+(n::Real, f::SisoRational) = t + n

-(f1::SisoRational, f2::SisoRational) = SisoRational(f1.num*f2.den - f2.num*f1.den, f1.den*f2.den)
-(n::Real, f::SisoRational) = SisoRational(n*f.den - f.num, f.den)
-(f::SisoRational, n::Real) = +(f, -n)
#.-(f::SisoRational, n::Real) = -(t, n)
#.-(n::Real, f::SisoRational) = -(n, t)

-(f::SisoRational) = SisoRational(-f.num, f.den)

*(f1::SisoRational, f2::SisoRational) = SisoRational(f1.num*f2.num, f1.den*f2.den)
*(f::SisoRational, n::Real) = SisoRational(f.num*n, f.den)
*(n::Real, f::SisoRational) = *(f, n)
#.*(f1::SisoRational, f2::SisoRational) = *(f1, f2)
#.*(f::SisoRational, n::Real) = *(f, n)
#.*(n::Real, f::SisoRational) = *(f, n)

/(n::Real, f::SisoRational) = SisoRational(n*f.den, f.num)
/(f::SisoRational, n::Real) = f*(1/n)
/(f1::SisoRational, f2::SisoRational) = f1*(1/f2)
