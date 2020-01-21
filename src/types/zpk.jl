""" `zpk(gain[, Ts]), zpk(num, den, k[, Ts]), zpk(sys)`

Create transfer function on zero pole gain form. The numerator and denominator are represented by their poles and zeros.

`sys::TransferFunction{SisoZpk{T,TR}} = k*numerator/denominator`
where `T` is the type of `k` and `TR` the type of the zeros/poles, usually Float64 and Complex{Float64}.

`num`: the roots of the numerator polynomial. Either scalar or vector to create SISO systems
or an array of vectors to create MIMO system.

`den`: the roots of the denominator polynomial. Either vector to create SISO systems
or an array of vectors to create MIMO system.

`k`: The gain of the system. Obs, this is not the same as `dcgain`.

`Ts`: Sample time if discrete time system.

Other uses:

`zpk(sys)`: Convert `sys` to `zpk` form.

`zpk("s")`: Create the transferfunction `s`.

"""
function zpk(z::VecOrMat{<:AbstractVector{TZ}}, p::VecOrMat{<:AbstractVector{TP}}, k::VecOrMat{T0}, Ts::SampleT) where {SampleT<:AbstractSampleTime, T0<:Number, TZ<:Number, TP<:Number}
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
    return TransferFunction{SampleT,SisoZpk{T,TR}}(matrix, Ts)
end
function zpk(z::AbstractVector{TZ}, p::AbstractVector{TP}, k::T, Ts::SampleT) where {SampleT<:AbstractSampleTime, T<:Number, TZ<:Number, TP<:Number}
    return zpk(fill(z, 1, 1), fill(p, 1, 1), fill(k, 1, 1), Ts)
end

function zpk(z::AbstractVector, p::AbstractVector, k::T, Ts::SampleT) where {SampleT<:AbstractSampleTime, T<:Number} # To be able to send in empty vectors [] of type Any
    if eltype(z) == Any && eltype(p) == Any
        @assert z == []
        @assert p == []
        return zpk(T[], T[], k, Ts)
    elseif eltype(z) == Any
        @assert z == []
        TR = eltype(p)
        return zpk(TR[], p, k, Ts)
    elseif eltype(p) == Any
        @assert p == []
        TR = eltype(z)
        return zpk(z, TR[], k, Ts)
    else
        error("Non numeric vectors must be empty.")
    end
end

function zpk(gain::Matrix{T}, Ts::SampleT; kwargs...) where {SampleT<:AbstractSampleTime, T <: Number}
    TR = promote_type(ComplexF64,T)
    ny, nu = size(gain)
    matrix = [SisoZpk{T, TR}(TR[],TR[], gain[o, i]) for o=1:ny, i=1:nu]
    return TransferFunction{SampleT,SisoZpk{T,TR}}(matrix, Ts)
end

zpk(k::Real, Ts::AbstractSampleTime) = zpk(eltype(k)[], eltype(k)[], k, Ts)


zpk(sys::StateSpace) = zpk(tf(sys)) # FIXME: probably better with direct conversion

function zpk(G::TransferFunction{SampleT,S}) where {SampleT<:AbstractSampleTime,T0, S<:SisoTf{T0}}
    T = typeof(one(T0)/one(T0))
    convert(TransferFunction{SampleT,SisoZpk{T, complex(T)}}, G)
end

zpk(var::AbstractString) = zpk(tf(var))
zpk(var::AbstractString, Ts::Real) = zpk(tf(var, Ts))

# Catch all 3(4) argument versions
zpk(z, p, k, Ts::Number) = zpk(z, p, k, Discrete(Ts))
zpk(z, p, k) = zpk(z, p, k, Continuous())
# Catch all 1(2) argument versions
zpk(gain, Ts::Number; kwargs...) where {T <: Number} = zpk(gain, Discrete(Ts))
zpk(gain; kwargs...) where {T <: Number} = zpk(gain, Continuous())
