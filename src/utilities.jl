# Misfit functions

## Predicates ##

@doc """`isstable(sys)`

Returns `true` if `sys` is stable, else returns `false`.""" ->
function isstable(sys::LTISystem)
    if sys.Ts == 0
        if any(real(pole(sys)).>0)
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
>= order(num))""" ->
function isproper(t::TransferFunction)
    for s in t.matrix
        if length(s.num) > length(s.den)
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

# Ensures the metadata for an LTISystem is valid
function validate_names(kwargs, key, n)
    names = get(kwargs, key, "")
    if names == ""
        return UTF8String[names for i = 1:n]
    elseif isa(names, Vector) && eltype(names) <: UTF8String
        return names
    elseif isa(names, UTF8String)
        return UTF8String[names * "$i" for i = 1:n]
    else
        error("$key must be of type `UTF8String` or Vector{UTF8String}")
    end
end

# Format the metadata for printing
function format_names(names::Vector{UTF8String}, default::String, unknown::String)
    n = size(names, 1)
    if all(names .== "")
        return UTF8String[default * string(i) for i=1:n]
    else
        for (i, n) in enumerate(names)
            names[i] = (n == "") ? unknown : n
        end
        return names
    end
end

# Creates an Array of size (ny,nx) with Vector elements from args
function array(T::Type, ny::Int,nx::Int, args::AbstractArray...)
  if ny*nx != length(args)
    error("Number of vectors must fit dimensions")
  end
  array = Array{Array{T,1},2}(ny,nx)

  for (i,v) in enumerate(args)
    row = floor(Int,(i-1)./nx)+1
    col = mod(i-1,nx)+1
    array[row,col] = v
  end
  array
end

                
array(ny::Int,nx::Int, args::AbstractArray...) = array(Float64, ny, nx, args...)
                  
