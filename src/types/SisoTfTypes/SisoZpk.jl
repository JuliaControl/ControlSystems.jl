# FIXME: Add some condition to guarantee that: TR <: Union(T, complex(T)) !

# T the numeric type of the transfer function
# TR the type of the roots
struct SisoZpk{T,TR<:Number} <: SisoTf{T}
    z::Vector{TR}
    p::Vector{TR}
    k::T
    function SisoZpk{T,TR}(z::Vector{TR}, p::Vector{TR}, k::T) where {T<:Number, TR<:Number}
        if k == zero(T)
            p = TR[]
            z = TR[]
        end
        if TR <: Complex && T <: Real
            z, p = copy(z), copy(p)
            @assert pairup_conjugates!(z) "zpk model should be real-valued, but zeros do not come in conjugate pairs."
            @assert pairup_conjugates!(p) "zpk model should be real-valued, but poles do not come in conjugate pairs."
        end
        new{T,TR}(z, p, k)
    end
end
function SisoZpk{T,TR}(z::Vector, p::Vector, k::Number) where {T<:Number, TR<:Number}
    SisoZpk{T,TR}(Vector{TR}(z), Vector{TR}(p), T(k))
end
function SisoZpk{T}(z::Vector, p::Vector, k::Number) where T
    TR = complex(T)
    SisoZpk{T,TR}(Vector{TR}(z), Vector{TR}(p), T(k))
end
function SisoZpk(z::AbstractVector{TZ}, p::AbstractVector{TP}, k::T) where {T<:Number, TZ<:Number, TP<:Number} # NOTE: is this constructor really needed?
    TR = promote_type(TZ,TP)
    SisoZpk{T,TR}(Vector{TR}(z), Vector{TR}(p), k)
end




Base.zero(::Type{SisoZpk{T}}) where T = SisoZpk{T}(T[], T[], zero(T))
Base.one(::Type{SisoZpk{T}}) where T = SisoZpk{T}(T[], T[], one(T))

Base.zero(::Type{SisoZpk{T,TR}}) where {T, TR} = SisoZpk{T,TR}(TR[], TR[], zero(T))
Base.one(::Type{SisoZpk{T,TR}}) where {T, TR} = SisoZpk{T,TR}(TR[], TR[], one(T))

Base.one(f::SisoZpk) = one(typeof(f))
Base.zero(f::SisoZpk) = zero(typeof(f))


# tzero is not meaningful for transfer function element? But both zero and zeros are taken...
tzero(f::SisoZpk) = f.z # Do minreal first?,
pole(f::SisoZpk) = f.p # Do minreal first?

numpoly(f::SisoZpk{<:Real}) = f.k*prod(roots2real_poly_factors(f.z))
denpoly(f::SisoZpk{<:Real}) = prod(roots2real_poly_factors(f.p))

numpoly(f::SisoZpk) = f.k*prod(roots2poly_factors(f.z))
denpoly(f::SisoZpk) = prod(roots2poly_factors(f.p))

numvec(f::SisoZpk) = reverse(coeffs(numpoly(f))) # FIXME: reverse?!
denvec(f::SisoZpk) = reverse(coeffs(denpoly(f))) # FIXME: reverse?!

isproper(f::SisoZpk) = (length(f.z) <= length(f.p))

# Will not be able to create complex polynomials without type instability


function minreal(sys::SisoZpk{T,TR}, eps::Real) where {T, TR}
    if length(sys.p) == 0
        return sys
    end

    newZ = copy(sys.z)
    newP = Vector{TR}()

    pidx = 1
    p = sys.p[pidx]
    while true
        if isempty(newZ)
            push!(newP, p)
        else
            distance, zidx = findmin(abs.(p .- newZ))

            if distance < eps
                if imag(p) == 0 && imag(newZ[zidx]) != 0
                    newZ[zidx+1] = real(newZ[zidx+1])
                end
                if imag(newZ[zidx]) == 0 && imag(p) != 0
                    pidx += 1
                    p = real(sys.p[pidx])
                    deleteat!(newZ, zidx)
                    continue
                end
                deleteat!(newZ, zidx)
            else
                push!(newP, p)
            end
        end

        pidx += 1
        if pidx > length(sys.p)
            break
        end
        p = sys.p[pidx]
    end
    SisoZpk{T, TR}(newZ, newP, sys.k)
end

""" Reorder the vector x of complex numbers so that complex conjugates come after each other, 
    with the one with positive imaginary part first. Returns true if the conjugates can be 
    paired and otherwise false."""
function pairup_conjugates!(x::AbstractVector)
    i = 0
    while i < length(x)
        i += 1
        imag(x[i]) == 0 && continue

        # Attempt to find a matching conjugate to x[i]
        j = findnext(==(conj(x[i])), x, i+1)
        j === nothing && return false

        tmp = x[j]
        x[j] = x[i+1]
        # Make sure that the complex number with positive imaginary part comes first
        if imag(x[i]) > 0
            x[i+1] = tmp
        else
            x[i+1] = x[i]
            x[i] = tmp
        end
        i += 1 # Since it is a pair and the conjugate was found
    end
    return true
end

function evalfr(f::SisoZpk{T1,TR}, s::T2) where {T1<:Number, TR<:Number, T2<:Number}
    T0 = promote_type(T2, TR)
    T = promote_type(T1, Base.promote_op(/, T0, T0))
    den = mapreduce(a -> s-a, *, f.p, init = one(T0)) # Default to one

    if den == zero(T0)
        return convert(T, Inf)
    else
        num = mapreduce(b -> s-b, *, f.z, init=one(T0))
        return f.k*num/den
    end
end


function poly_factors2string(poly_factors::AbstractArray{<:Polynomial{T}}, var) where T
    if length(poly_factors) == 0
        str = sprint(printpolyfun(var), Polynomial(one(T)))
    elseif length(poly_factors) == 1
        str = sprint(printpolyfun(var), poly_factors[1])
    else
        str = prod(z -> "("*sprint(printpolyfun(var), z)*")", poly_factors)
    end
end

""" Heurisitc function that tries to add parentheses when printing the coeffcient
    for systems on zpk form. Should at least handle the following types
    Measurment, Dual, Sym. """
function _printcoefficient(nbr::Number)
    # Print type information as in 1.0f0 for Float32
    # showcompact might be better, but is not consistent with polynomials
    nbr_string = sprint(show,nbr)
    if occursin(" + ", nbr_string) || occursin(" - ", nbr_string) || occursin(" ± ", nbr_string)
        return "(" * nbr_string * ")" # Add parens
    else
        return nbr_string
    end
end

function print_siso(io::IO, t::SisoZpk, var=:s)
    if typeof(t.k) <: Real
        numstr = poly_factors2string(roots2real_poly_factors(t.z), var)
        denstr = poly_factors2string(roots2real_poly_factors(t.p), var)
    else
        numstr = poly_factors2string(roots2poly_factors(t.z), var)
        denstr = poly_factors2string(roots2poly_factors(t.p), var)
    end
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

    gainstr = _printcoefficient(t.k)
    #Add spaces to account for gain string
    numstr = " "^(length(gainstr))*numstr
    denstr = " "^(length(gainstr))*denstr
    println(io, numstr)
    println(io, gainstr*repeat("-", dashcount))
    println(io, denstr)
end





==(f1::SisoZpk, f2::SisoZpk) = (f1-f2).k == 0
function isapprox(f1::SisoZpk, f2::SisoZpk; rtol::Real=sqrt(eps()), atol::Real=sqrt(eps()))
    fdiff = f1 - f2
    isapprox(fdiff.k, 0, atol=atol, rtol=rtol)
end

function +(f1::SisoZpk{T1,TR1}, f2::SisoZpk{T2,TR2}) where {T1<:Number,T2<:Number,TR1<:Number,TR2<:Number}
    numPoly = numpoly(f1)*denpoly(f2) + numpoly(f2)*denpoly(f1)

    TRtmp = promote_type(TR1, TR2)
    # Calculating roots can make integers floats
    TRnew = Base.promote_op(/, TRtmp, TRtmp)
    z = convert(Vector{TRnew}, roots(numPoly))
    #TODO gains could remain integers, but numerical precision inhibits this
    Ttmp = promote_type(T1,T2)
    Tnew = Base.promote_op(/, Ttmp, Ttmp)
    if length(numPoly) > 0
        k = convert(Tnew, numPoly[end])
        p = convert(Vector{TRnew}, [f1.p;f2.p])
    else
        k = zero(Tnew)
        p = TRnew[]
    end

    # FIXME:
    # Threshold for pole-zero cancellation should depend on the roots of the system
    # Note the difference between continuous and discrete-time systems...
    minreal(SisoZpk(z::Vector{TRnew},p::Vector{TRnew},k), sqrt(eps())) # type assert required or inference
end


+(f::SisoZpk, n::T) where {T<:Number} = f + SisoZpk{T,T}(T[],T[],n)
+(n::Number, f::SisoZpk) = f + n
#.+(t::SisoZpk, n::Number) = t + n
#.+(n::Number, t::SisoZpk) = t + n

-(t1::SisoZpk, t2::SisoZpk) = +(t1,-t2)
-(n::T, f::SisoZpk) where {T<:Number} = SisoZpk{T,T}(T[],T[],n) - f
-(t::SisoZpk, n::Number) = +(t, -n)
#.-(t::SisoZpk, n::Number) = -(t, n)
#.-(n::Number, t::SisoZpk) = -(n, t)


-(f::SisoZpk) = SisoZpk(f.z, f.p, -f.k)


*(f1::SisoZpk, f2::SisoZpk) = SisoZpk([f1.z;f2.z], [f1.p;f2.p], f1.k*f2.k)

*(f::SisoZpk, n::Number) = SisoZpk(f.z, f.p, f.k*n)
*(n::Number, f::SisoZpk) = *(f, n)
#.*(f1::SisoZpk, f2::SisoZpk) = *(f1, f2)
#.*(t::SisoZpk, n::Number) = *(t, n)
#.*(n::Number, t::SisoZpk) = *(t, n)

/(n::Number, f::SisoZpk) = SisoZpk(f.p, f.z, n/f.k)
/(f::SisoZpk, n::Number) = SisoZpk(f.z, f.p, f.k/n)
/(f1::SisoZpk, f2::SisoZpk) = f1*(1/f2)
