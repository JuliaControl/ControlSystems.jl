# T the numeric type of the transfer function
# TR the type of the roots

# NOTE: Real is not a subtype of complex...

# Add some condition like: TR <: Union(T, complex(T)) ?

struct SisoZpk{T,TR<:Number} <: SisoTf{T}
    z::Vector{TR}
    p::Vector{TR}
    k::T
    function SisoZpk{T,TR}(z::Vector{TR}, p::Vector{TR}, k::T) where  {T<:Number, TR<:Number}
        if k == zero(T)
            p = TR[]
            z = TR[]
        end
        new{T,TR}(z, p, k)
    end
end
function SisoZpk{T,TR}(z::Vector, p::Vector, k::Number) where  {T<:Number, TR<:Number}
    SisoZpk{T,TR}(Vector{TR}(z), Vector{TR}(p), T(k))
end
function SisoZpk{T}(z::Vector, p::Vector, k::Number) where T
    TR = complex(T)
    SisoZpk{T,TR}(Vector{TR}(z), Vector{TR}(p), T(k))
end
function SisoZpk(z::AbstractArray{TZ}, p::AbstractArray{TP}, k::T) where {T<:Number, TZ<:Number, TP<:Number} # NOTE: is this constructor really needed?
      TR = promote_type(TZ,TP,T)
      # Should check if roots matches up

      if TR <: Complex && T <: Real
          check_real # Only throws an error
      end
      SisoZpk{T,TR}(Vector{TR}(z), Vector{TR}(p), k)
end

Base.zero(::Type{SisoZpk{T}}) where T = SisoZpk{T}(T[], T[], zero(T))
Base.one(::Type{SisoZpk{T}}) where T = SisoZpk{T}(T[], T[], one(T))

Base.zero(::Type{SisoZpk{T,TR}}) where {TR, T<:TR} = SisoZpk{T,TR}(TR[], TR[], zero(T))
Base.one(::Type{SisoZpk{T,TR}}) where {TR, T<:TR} = SisoZpk{T,TR}(TR[], TR[], one(T))

Base.one(f::SisoZpk) = one(typeof(f))
Base.zero(f::SisoZpk) = zero(typeof(f))


# tzero is not meaningful for transfer function element? But both zero and zeros are taken...
tzero(sys::SisoZpk) = f.z # Do minreal first?,
pole(sys::SisoZpk) = f.p # Do minreal first?

numpoly(f::SisoZpk) = prod(roots2real_poly_factors(f.z))
denpoly(f::SisoZpk) = prod(roots2real_poly_factors(f.p))

num(f::SisoZpk) = reverse(coeffs(numpoly(f))) # FIXME: reverse?!
den(f::SisoZpk) = reverse(coeffs(denpoly(f))) # FIXME: reverse?!


# Will not be able to create complex polynomials without type instability


# If TR is Complex and T is Real, check that every pole is matched to its conjugate
function check_real(r_vec::AbstractVector{Complex})
    for k=1:length(r_vec)
        if isreal(r_vec[k])
            continue
        elseif k == length(r_vec)
            error("No match found")
        elseif conj(r_vec[k]) == r_vec[k+1]
            k=k+1
            continue
        else #
            # findfirst(conj(r_vec[k]), r_vec[k+1:end])
            # if no element found, throw error
            # r_vec[k+1] = conjuagte
        end
    end
end



function evalfr(sys::SisoZpk, s::Number)
    S = promote_type(typeof(s), Float64)
    den = prod(s-a for a in sys.p)
    if den == zero(S)
        return convert(S, Inf)
    else
        num = prod(s-b for b in sys.z)
        num, typeof(num), zp2polys(sys.z)
        return sys.k*num/den
    end
end

function poly_factors2string(poly_factors)
    if length(poly_factors) == 0
        str = "1.0"
    elseif length(poly_factors) == 1
        str = sprint(printpoly, poly_factors[1])
    else
        str = reduce(*,"",["("*sprint(printpoly, z)*")" for z in poly_factors])
    end
end

# TODO: add print function for complex coefficient polynomial factors
function print_siso(io::IO, t::SisoZpk, var=:s)
    numstr = poly_factors2string(roots2real_poly_factors(t.z))
    denstr = poly_factors2string(roots2real_poly_factors(t.p))

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

    gainstr = string(t.k)
    #Add spaces to account for gain string
    numstr = " "^(length(gainstr))*numstr
    denstr = " "^(length(gainstr))*denstr
    println(io, numstr)
    println(io, gainstr*repeat("-", dashcount))
    println(io, denstr)
end





==(t1::SisoZpk, t2::SisoZpk) = (t1-t2).k == 0.0
function isapprox(t1::SisoZpk, t2::SisoZpk; rtol::Real=sqrt(eps()), atol::Real=sqrt(eps()))
    tdiff = t1 - t2
    isapprox(tdiff.k, 0, atol=atol, rtol=rtol)
end

function +(t1::SisoZpk, t2::SisoZpk)
  numPoly = t1.k*prod(zp2polys(t1.z))*prod(zp2polys(t2.p))+t2.k*prod(zp2polys(t2.z))*prod(zp2polys(t1.p))
  z = roots(numPoly)
  if length(numPoly) > 0
      k = numPoly[1]
      p = [t1.p;t2.p]
  else
      k = 0
      p = []
  end
  SisoZpk(z,p,k)
end

+(t::SisoZpk, n::Real) = t + SisoZpk([],[],n)
+(n::Real, t::SisoZpk) = SisoZpk([],[],n) + t
#.+(t::SisoZpk, n::Real) = t + n
#.+(n::Real, t::SisoZpk) = t + n

-(t1::SisoZpk, t2::SisoZpk) = +(t1,-t2)
-(n::Real, t::SisoZpk) = SisoZpk([],[],n) - t
-(t::SisoZpk, n::Real) = +(t, -n)
#.-(t::SisoZpk, n::Real) = -(t, n)
#.-(n::Real, t::SisoZpk) = -(n, t)

-(t::SisoZpk) = SisoZpk(t.z, t.p, -t.k)

*(t1::SisoZpk, t2::SisoZpk) = SisoZpk([t1.z;t2.z], [t1.p;t2.p], t1.k*t2.k)
*(t::SisoZpk, n::Real) = SisoZpk(t.z, t.p, t.k*n)
*(n::Real, t::SisoZpk) = *(t, n)
#.*(t1::SisoZpk, t2::SisoZpk) = *(t1, t2)
#.*(t::SisoZpk, n::Real) = *(t, n)
#.*(n::Real, t::SisoZpk) = *(t, n)

/(n::Real, t::SisoZpk) = SisoZpk(t.p, t.z, n/t.k)
/(t::SisoZpk, n::Real) = SisoZpk(t.z, t.p, t.k/n)
/(t1::SisoZpk, t2::SisoZpk) = t1*(1/t2)
