""" `zpk(gain, Ts=0), zpk(num, den, k, Ts=0), zpk(sys)`

Create transfer function on zero pole gain form. The numerator and denominator are represented by their poles and zeros.

`sys::TransferFunction{SisoZpk{T,TR}} = k*numerator/denominator`
where `T` is the type of `k` and `TR` the type of the zeros/poles, usually Float64 and Complex{Float64}.

`num`: the roots of the numerator polynomial. Either scalar or vector to create SISO systems
or an array of vectors to create MIMO system.

`den`: the roots of the denominator polynomial. Either vector to create SISO systems
or an array of vectors to create MIMO system.

`k`: The gain of the system. Obs, this is not the same as `dcgain`.

`Ts`: Sample time or `0` for continuous system.

Other uses:

`zpk(sys)`: Convert `sys` to `zpk` form.

`zpk("s")`: Create the transferfunction `s`.

"""
function zpk(z::VecOrMat{<:AbstractVector{TZ}}, p::VecOrMat{<:AbstractVector{TP}}, k::VecOrMat{T0}, Ts::Real=0.0) where {T0<:Number, TZ<:Number, TP<:Number}
    # Validate input and output dimensions match
    if !(size(z) == size(p) == size(k))
        error("Dimensions for z, p, and k must match")
    end

    TR = promote_type(T0,TZ,TP)
    T = promote_type(T0, real(TR)) # Ensuring this avoids many problems, e.g., with SisoZpk{Int64,ComplexF64}

    matrix = Array{SisoZpk{T, TR}}(undef, size(z,1), size(z,2)) # TODO: make this nicer
    for o=1:size(z,1)
        for i=1:size(z,2)
            matrix[o, i] = SisoZpk{T,TR}(z[o, i], p[o, i], k[o, i])
        end
    end
    return TransferFunction{SisoZpk{T,TR}}(matrix, Float64(Ts))
end
function zpk(z::AbstractVector{TZ}, p::AbstractVector{TP}, k::T, Ts::Real=0.0) where {T<:Number, TZ<:Number, TP<:Number}
    return zpk(fill(z, 1, 1), fill(p, 1, 1), fill(k, 1, 1), Ts)
end

function zpk(gain::Matrix{T}, Ts::Real=0; kwargs...) where {T <: Number}
    TR = promote_type(ComplexF64,T)
    ny, nu = size(gain)
    matrix = [SisoZpk{T, TR}(TR[],TR[], gain[o, i]) for o=1:ny, i=1:nu]
    return TransferFunction{SisoZpk{T,TR}}(matrix, Ts)
end

function zpk(z::AbstractVector, p::AbstractVector, k::T) where {T<:Number} # To be able to send in empty vectors [] of type Any
    if eltype(z) == Any && eltype(p) == Any
        @assert z == []
        @assert p == []
        return zpk(T[], T[], k)
    elseif eltype(z) == Any
        @assert z == []
        TR = eltype(p)
        return zpk(TR[], p, k)
    elseif eltype(p) == Any
        @assert p == []
        TR = eltype(z)
        return zpk(z, TR[], k)
    else
        error("Non numeric vectors must be empty.")
    end
end

zpk(k::Real, Ts::Real=0) = zpk(eltype(k)[], eltype(k)[], k, Ts)

zpk(sys::StateSpace) = zpk(tf(sys)) # FIXME: probably better with direct conversion

function zpk(G::TransferFunction{S}) where {T0, S<:SisoTf{T0}}
    T = typeof(one(T0)/one(T0))
    convert(TransferFunction{SisoZpk{T, complex(T)}}, G)
end

zpk(var::AbstractString) = zpk(tf(var))
zpk(var::AbstractString, Ts::Real) = zpk(tf(var, Ts))
