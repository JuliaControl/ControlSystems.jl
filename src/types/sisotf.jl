## User should just use TransferFunction
immutable SisoTf
    num::Poly{Float64}
    den::Poly{Float64}
    function SisoTf(num::Poly{Float64}, den::Poly{Float64})
        if all(den == zero(den))
            error("Zero denominator")
        elseif all(num == zero(num))
            # The numerator is zero, make the denominator 1
            den = one(den)
        end
        new(num, den)
    end
end
SisoTf(num::Vector, den::Vector) = SisoTf(Poly(map(Float64,num)), Poly(map(Float64,den)))

function print_sisotf(io::IO, t::SisoTf, var=:s)
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

Base.promote_rule{T<:Real}(::Type{SisoTf}, ::Type{T}) = SisoTf
Base.convert(::Type{SisoTf}, b::Real) = SisoTf([b], [1])

Base.zero(::Type{SisoTf}) = SisoTf(zero(Poly{Float64}), one(Poly{Float64}))
Base.zero(t::SisoTf) = Base.zero(SisoTf)

Base.length(t::SisoTf) = max(length(t.num), length(t.den))

function Base.num(t::SisoTf)
    lt = length(t)
    n = zeros(lt)
    n[(lt - length(t.num) + 1):end] = t.num[:]
    return n
end

function Base.den(t::SisoTf)
    lt = length(t)
    d = zeros(lt)
    d[(lt - length(t.den) + 1):end] = t.den[:]
    return d
end

==(t1::SisoTf, t2::SisoTf) = (t1.num == t2.num && t1.den == t2.den) 

+(t1::SisoTf, t2::SisoTf) = SisoTf(t1.num*t2.den + t2.num*t1.den, t1.den*t2.den)
+(t::SisoTf, n::Real) = SisoTf(t.num + n*t.den, t.den)
+(n::Real, t::SisoTf) = t + n
.+(t::SisoTf, n::Real) = t + n
.+(n::Real, t::SisoTf) = t + n

-(t1::SisoTf, t2::SisoTf) = SisoTf(t1.num*t2.den - t2.num*t1.den, t1.den*t2.den)
-(n::Real, t::SisoTf) = SisoTf(n*t.den - t.num, t.den)
-(t::SisoTf, n::Real) = +(t, -n)
.-(t::SisoTf, n::Real) = -(t, n)
.-(n::Real, t::SisoTf) = -(n, t)

-(t::SisoTf) = SisoTf(-t.num, t.den)

*(t1::SisoTf, t2::SisoTf) = SisoTf(t1.num*t2.num, t1.den*t2.den)
*(t::SisoTf, n::Real) = SisoTf(t.num*n, t.den)
*(n::Real, t::SisoTf) = *(t, n)
.*(t1::SisoTf, t2::SisoTf) = *(t1, t2)
.*(t::SisoTf, n::Real) = *(t, n)
.*(n::Real, t::SisoTf) = *(t, n)

/(n::Real, t::SisoTf) = SisoTf(n*t.den, t.num)
/(t::SisoTf, n::Real) = t*(1/n)
/(t1::SisoTf, t2::SisoTf) = t1*(1/t2)
