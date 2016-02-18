## User should just use TransferFunction
immutable SisoRational <: SisoTf
    num::Poly{Float64}
    den::Poly{Float64}
    function SisoRational(num::Poly{Float64}, den::Poly{Float64})
        if all(den == zero(den))
            error("Zero denominator")
        elseif all(num == zero(num))
            # The numerator is zero, make the denominator 1
            den = one(den)
        end
        new(num, den)
    end
end
SisoRational(num::Vector, den::Vector) = SisoRational(Poly(map(Float64,num)), Poly(map(Float64,den)))

function print_sisorational(io::IO, t::SisoRational, var=:s)
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

Base.promote_rule{T<:Real}(::Type{SisoRational}, ::Type{T}) = SisoRational
Base.convert(::Type{SisoRational}, b::Real) = SisoRational([b], [1])

Base.zero(::Type{SisoRational}) = SisoRational(zero(Poly{Float64}), one(Poly{Float64}))
Base.zero(t::SisoRational) = Base.zero(SisoRational)

Base.length(t::SisoRational) = max(length(t.num), length(t.den))

function Base.num(t::SisoRational)
    lt = length(t)
    n = zeros(lt)
    n[(lt - length(t.num) + 1):end] = t.num[:]
    return n
end

function Base.den(t::SisoRational)
    lt = length(t)
    d = zeros(lt)
    d[(lt - length(t.den) + 1):end] = t.den[:]
    return d
end

==(t1::SisoRational, t2::SisoRational) = (t1.num == t2.num && t1.den == t2.den)

+(t1::SisoRational, t2::SisoRational) = SisoRational(t1.num*t2.den + t2.num*t1.den, t1.den*t2.den)
+(t::SisoRational, n::Real) = SisoRational(t.num + n*t.den, t.den)
+(n::Real, t::SisoRational) = t + n
.+(t::SisoRational, n::Real) = t + n
.+(n::Real, t::SisoRational) = t + n

-(t1::SisoRational, t2::SisoRational) = SisoRational(t1.num*t2.den - t2.num*t1.den, t1.den*t2.den)
-(n::Real, t::SisoRational) = SisoRational(n*t.den - t.num, t.den)
-(t::SisoRational, n::Real) = +(t, -n)
.-(t::SisoRational, n::Real) = -(t, n)
.-(n::Real, t::SisoRational) = -(n, t)

-(t::SisoRational) = SisoRational(-t.num, t.den)

*(t1::SisoRational, t2::SisoRational) = SisoRational(t1.num*t2.num, t1.den*t2.den)
*(t::SisoRational, n::Real) = SisoRational(t.num*n, t.den)
*(n::Real, t::SisoRational) = *(t, n)
.*(t1::SisoRational, t2::SisoRational) = *(t1, t2)
.*(t::SisoRational, n::Real) = *(t, n)
.*(n::Real, t::SisoRational) = *(t, n)

/(n::Real, t::SisoRational) = SisoRational(n*t.den, t.num)
/(t::SisoRational, n::Real) = t*(1/n)
/(t1::SisoRational, t2::SisoRational) = t1*(1/t2)
