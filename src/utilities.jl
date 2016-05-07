# Misfit functions

## Predicates ##

@doc """`isstable(sys)`

Returns `true` if `sys` is stable, else returns `false`.""" ->
function isstable(sys::LTISystem)
    if sys.Ts == 0
        if any(real(pole(sys)).>=0)
            return false
        end
    else
        if any(abs(pole(sys)).>=1)
            return false
        end
    end
    return true
end

@doc """`iscontinuous(sys)`

Returns `true` if `sys` is continuous, else returns `false`.""" ->
function iscontinuous(sys::LTISystem)
    return sys.Ts == 0
end

@doc """`issiso(sys)`

Returns `true` if `sys` is SISO, else returns `false`.""" ->
function issiso(sys::LTISystem)
    return sys.ny == 1 && sys.nu == 1
end

@doc """`isproper(tf)`

Returns `true` if the `TransferFunction` is proper. This means that order(den)
\>= order(num))""" ->
function isproper(t::TransferFunction)
    for s in t.matrix
        if length(numpoly(s)) > length(denpoly(s))
            return false
        end
    end
    return true
end

## Helpers (not public) ##

# Convert the argument to a Matrix{Float64}
function float64mat(A::Vector)
    A = reshape(A, size(A, 1), 1)
    return float64mat(A)
end
float64mat(A::Matrix) = map(Float64,A)
float64mat(A::Matrix{Float64}) = A
float64mat(A::Real) = float64mat([A])

# Ensures the metadata for an LTISystem is valid
function validate_names(kwargs, key, n)
    names = get(kwargs, key, "")
    if names == ""
        return fill(UTF8String(""), n)
    elseif isa(names, Vector) && eltype(names) <: AbstractString
        return UTF8String[names[i] for i = 1:n]
    elseif isa(names, AbstractString)
        return UTF8String[names * "$i" for i = 1:n]
    else
        error("$key must be of type `AbstractString` or Vector{AbstractString}")
    end
end

# Format the metadata for printing
function format_names(names::Vector{UTF8String}, default::AbstractString, unknown::AbstractString)
    n = size(names, 1)
    if all(names .== "")
        return UTF8String[default * string(i) for i=1:n]
    else
        for (i, n) in enumerate(names)
            names[i] = UTF8String((n == "") ? unknown : n)
        end
        return names
    end
end

function unwrap!(M::Array, dim=1)
    alldims(i) = [ n==dim ? i : (1:size(M,n)) for n in 1:ndims(M) ]
    for i = 2:size(M, dim)
        #This is a copy of slicedim from the JuliaLang but enables us to write to it
        #The code (with dim=1) is equivalent to
        # d = M[i,:,:,...,:] - M[i-1,:,...,:]
        # M[i,:,:,...,:] -= floor((d+π) / (2π)) * 2π
        d = M[alldims(i)...] - M[alldims(i-1)...]
        M[alldims(i)...] -= floor((d+π) / 2π) * 2π
    end
    return M
end

#Collect will create a copy and collect the elements
unwrap(m::AbstractArray, args...) = unwrap!(collect(m), args...)
unwrap(x::Number) = x

poly2vec(p::Poly) = p.a[1:end]
