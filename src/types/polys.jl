# Taken and modified from https://github.com/vtjnash/Polynomial.jl
# Polynomial.jl was deprecated in favor of Polynomials.jl, which uses reversed
# indexing. The former was copied to this file, and modified/trimmed to just
# the bare functions required for TransferFunction support.

Base.eps{T}(::Type{T}) = zero(T)
Base.eps{F<:AbstractFloat}(x::Type{F}) = Base.eps(F)
Base.eps{T}(x::Type{Complex{T}}) = eps(T)

immutable Poly{T<:Number}
    a::Vector{T}
    nzfirst::Int #for effiencicy, track the first non-zero index
    function Poly(a::Vector{T})
        la = length(a)
        i = 0
        for i = 1:la
            if abs(a[i]) > 2*eps(T)  break  end
        end
        new(a, i)
    end
end

Poly{T<:Number}(a::Vector{T}) = Poly{T}(a)

Base.convert{T}(::Type{Poly{T}}, p::Poly) = Poly(convert(Vector{T}, p.a))
Base.promote_rule{T, S}(::Type{Poly{T}}, ::Type{Poly{S}}) = Poly{promote_type(T, S)}
Base.eltype{T}(::Poly{T}) = T

Base.length(p::Poly) = length(p.a) - p.nzfirst + 1
Base.endof(p::Poly) = length(p)
deg(p::Poly) = length(p) - 1

Base.getindex(p::Poly, i) = p.a[i - 1 + p.nzfirst]
Base.getindex(p::Poly, c::Colon) = p.a[p.nzfirst:end]
Base.setindex!(p::Poly, v, i) = (p.a[i - 1 + p.nzfirst] = v)

Base.copy(p::Poly) = Poly(copy(p.a[p.nzfirst:end]))

Base.zero{T}(p::Poly{T}) = Poly([zero(T)])
Base.zero{T}(::Type{Poly{T}}) = Poly([zero(T)])
Base.one{T}(p::Poly{T}) = Poly([one(T)])
Base.one{T}(::Type{Poly{T}}) = Poly([one(T)])

function Base.show(io::IO, p::Poly)
    print(io,"Poly(")
    print(io,p)
    print(io,")")
end

Base.print(io::IO, p::Poly) = print_poly(io, p)

function print_poly{T}(io::IO, p::Poly{T}, var=:x)
    n = length(p)
    if n == 1
        print(io, p[1])
    else
        for j = 1:n
            pj = p[j]
            magpj = abs(pj)
            if magpj > 2*eps(T)
                if j == 1
                    real(pj) < 0 && print(io, "-")    #Prepend - if first and negative
                else
                    real(pj) < 0 ? print(io," - ") : print(io," + ")
                end
                #Print pj if pj is the last coefficient, or pj is not identically 1
                if j == n || abs(magpj - 1) > 2*eps(T)
                    print(io, magpj)
                end
                exp = n-j
                if exp > 0
                    print(io, var)
                    if exp > 1
                        print(io, '^', exp)
                    end
                end
            end
        end
    end
end

*(c::Number, p::Poly) = Poly(c * p.a[p.nzfirst:end])
*(p::Poly, c::Number) = Poly(c * p.a[p.nzfirst:end])
/(p::Poly, c::Number) = Poly(p.a[p.nzfirst:end] / c)
./(p::Poly, c::Number) = /(p, c)
-(p::Poly) = Poly(-p.a[p.nzfirst:end])

-(p::Poly, c::Number) = +(p, -c)
+(c::Number, p::Poly) = +(p, c)
function +(p::Poly, c::Number)
    if length(p) < 1
        return Poly([c,])
    else
        p2 = copy(p)
        p2.a[end] += c;
        return p2;
    end
end
function -(c::Number, p::Poly)
    if length(p) < 1
        return Poly([c,])
    else
        p2 = -p;
        p2.a[end] += c;
        return p2;
    end
end

function +{T,S}(p1::Poly{T}, p2::Poly{S})
    R = promote_type(T,S)
    n = length(p1)
    m = length(p2)
    if n > m
        a = Array(R, n)
        for i = 1:m
            a[n-m+i] = p1[n-m+i] + p2[i]
        end
        for i = 1:n-m
            a[i] = p1[i]
        end
    else
        a = Array(R, m)
        for i = 1:n
            a[m-n+i] = p1[i] + p2[m-n+i]
        end
        for i = 1:m-n
            a[i] = p2[i]
        end
    end
    Poly(a)
end

function -{T,S}(p1::Poly{T}, p2::Poly{S})
    R = promote_type(T,S)
    n = length(p1)
    m = length(p2)
    if n > m
        a = Array(R, n)
        for i = 1:m
            a[n-m+i] = p1[n-m+i] - p2[i]
        end
        for i = 1:n-m
            a[i] = p1[i]
        end
    else
        a = Array(R, m)
        for i = 1:n
            a[m-n+i] = p1[i] - p2[m-n+i]
        end
        for i = 1:m-n
            a[i] = -p2[i]
        end
    end
    Poly(a)
end

function *{T,S}(p1::Poly{T}, p2::Poly{S})
    R = promote_type(T,S)
    n = length(p1)
    m = length(p2)
    if n == 0 || m == 0
        return Poly(R[])
    end
    a = zeros(R, n+m-1)
    for i = 1:length(p1)
        for j = 1:length(p2)
            a[i+j-1] += p1[i] * p2[j]
        end
    end
    Poly(a)
end

function ==(p1::Poly, p2::Poly)
    if length(p1) != length(p2)
        return false
    else
        return p1.a[p1.nzfirst:end] == p2.a[p2.nzfirst:end]
    end
end

function isapprox(p1::Poly, p2::Poly; kwargs...)
    if length(p1) != length(p2)
        return false
    else
        return isapprox(p1.a[p1.nzfirst:end], p2.a[p2.nzfirst:end]; kwargs...)
    end
end

function polyval{T}(p::Poly{T}, x::Number)
    R = promote_type(T, typeof(x))
    lenp = length(p)
    if lenp == 0
        return zero(R)
    else
        y = convert(R, p[1])
        for i = 2:lenp
            y = p[i] + x.*y
        end
        return y
    end
end

Base.call{T}(p::Poly{T}, x::Number) = polyval(p, x)

# compute the roots of a polynomial
function roots{T}(p::Poly{T})
    R = promote_type(T, Float64)
    num_zeros = 0
    if length(p) == 0
        return zeros(R, 0)
    end
    while abs(p[end-num_zeros]) <= 2*eps(T)
        if num_zeros == length(p)-1
            return zeros(R, 0)
        end
        num_zeros += 1
    end
    n = length(p)-num_zeros-1
    if n < 1
        return zeros(R, length(p)-1)
    end
    companion = zeros(R, n, n)
    a0 = p[end-num_zeros]
    for i = 1:n-1
        companion[1,i] = -p[end-num_zeros-i] / a0
        companion[i+1,i] = 1;
    end
    companion[1,end] = -p[1] / a0
    D,V = eig(companion)
    r = zeros(eltype(D),length(p)-1)
    r[1:n] = 1./D
    return r
end
