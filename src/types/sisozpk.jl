RealOrComplex = Union{Real,Complex}

## User should just use TransferFunction
immutable SisoZpk
    z::Vector{Complex{Float64}}
    p::Vector{Complex{Float64}}
    k::Float64
    function SisoZpk(num::Vector{Complex{Float64}}, den::Vector{Complex{Float64}}, k::Float64)
        if k == zero(k)
            # The numerator is zero, make the denominator 1
            den = []
        end
        new(num, den, k)
    end
end

SisoZpk{T<:RealOrComplex,S<:RealOrComplex}(num::Vector{T}, den::Vector{S}, K::Real) = SisoZpk(float(num), float(den), float(K))


function zp2polys(vec)
    polys = Array{Poly{Float64},1}
    polesLeft = Set(1:length(vec))
    while length(polesLeft) > 0
        p = polys[pop!(polesLeft)]
        if abs(imag(p) < sqrt(eps()))
            push!(polys,Poly(float([1, -real(p)])))
        else
            polesTest = Complex128[i for i in s]
            val, i = findmin(abs(polesTest-conj(p)))
            val > sqrt(eps) && error("Could not find conjugate to pole")
            push!(polys,Poly(float([real(p)^2, 2*real(p)*imag(p), imag(p)^2])))
            pop!(polesLeft,polesTest[i])
        end
    end
end

function print_sisozpk(io::IO, t::SisoZpk, var=:s)
    zpolys = zp2polys(t.z)
    ppolys = zp2polys(t.p)
    # Convert the numerator and denominator to strings
    numstr = prod(["("*sprint(print_poly, z, var)*")" for z in zpolys])
    denstr = prod(["("*sprint(print_poly, p, var)*")" for p in ppolys])

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
    
    gainstr = @sprintf "%3.2f" t.k
    #Add spaces to account for gain string
    numstr = " "^(length(gainstr))*numstr
    denstr = " "^(length(gainstr))*denstr
    println(io, numstr)
    println(io, gainstr*repeat("-", dashcount))
    println(io, denstr)
end

Base.promote_rule{T<:Real}(::Type{SisoZpk}, ::Type{T}) = SisoZpk
Base.convert(::Type{SisoZpk}, b::Real) = SisoZpk([], [], b)

Base.zero(::Type{SisoZpk}) = SisoZpk([],[],0.0)
Base.zero(t::SisoZpk) = Base.zero(SisoTf)

Base.length(t::SisoZpk) = max(length(t.z), length(t.p))

# function Base.num(t::SisoZpk)
#     lt = length(t)
#     n = zeros(lt)
#     n[(lt - length(t.num) + 1):end] = t.num[:]
#     return n
# end
# 
# function Base.den(t::SisoTf)
#     lt = length(t)
#     d = zeros(lt)
#     d[(lt - length(t.den) + 1):end] = t.den[:]
#     return d
# end

==(t1::SisoZpk, t2::SisoZpk) = (t1.z == t2.z && t1.p == t2.p && t1.k == t2.k)

function +(t1::SisoZpk, t2::SisoZpk)
  #TODO Make sure gain is computed correctly
  num = poles(t1.k*prod(zp2polys(t1.num))*prod(zp2polys(t2.den))+t2.k*prod(zp2polys(t2.num))*prod(zp2polys(t1.den)))
  SisoZpk(num,[t1.den;t2.den],1)
end
+(t::SisoZpk, n::Real) = t + SisoZpk([],[],n)
+(n::Real, t::SisoZpk) = SisoZpk([],[],n) + t
.+(t::SisoZpk, n::Real) = t + n
.+(n::Real, t::SisoZpk) = t + n

-(t1::SisoZpk, t2::SisoZpk) = +(t1,-t2)
-(n::Real, t::SisoZpk) = SisoZpk([],[],n) - t
-(t::SisoZpk, n::Real) = +(t, -n)
.-(t::SisoZpk, n::Real) = -(t, n)
.-(n::Real, t::SisoZpk) = -(n, t)

-(t::SisoZpk) = SisoZpk(t.z, t.p, -z.k)

*(t1::SisoZpk, t2::SisoZpk) = SisoZpk([t1.z;t2.z], [t1.p;t2.z], t1.k*t2.k)
*(t::SisoZpk, n::Real) = SisoZpk(t.z, t.p, t.k*n)
*(n::Real, t::SisoZpk) = *(t, n)
.*(t1::SisoZpk, t2::SisoZpk) = *(t1, t2)
.*(t::SisoZpk, n::Real) = *(t, n)
.*(n::Real, t::SisoZpk) = *(t, n)

/(n::Real, t::SisoZpk) = SisoZpk(t.p, t.z, 1/t.k)
/(t::SisoZpk, n::Real) = SisoZpk(t.z, t.p, 1/t.k)
/(t1::SisoZpk, t2::SisoZpk) = t1*(1/t2)
