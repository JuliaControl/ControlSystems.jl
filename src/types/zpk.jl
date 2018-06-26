@doc """ `zpk(gain, Ts=0; kwargs...), zpk(num, den, k, Ts=0; kwargs...), zpk(sys)`

Create transfer function on zero pole gain form. The numerator and denominator are represented by their poles and zeros.

`sys = k*numerator/denominator`

`num`: the roots of the numerator polynomial. Either scalar or vector to create SISO systems
or an array of vectors to create MIMO system.

`den`: the roots of the denominator polynomial. Either vector to create SISO systems
or an array of vectors to create MIMO system.

`k`: The gain of the system. Obs, this is not the same as `dcgain`.

`Ts`: Sample time or `0` for continuous system.

Other uses:

`tf(sys)`: Convert `sys` to `tf` form.

`tf("s")`: Create the transferfunction `s`.""" ->
function zpk(z::VecOrMat{<:AbstractVector{TZ}}, p::VecOrMat{<:AbstractVector{TP}}, k::VecOrMat{T}, Ts::Real=0.0) where {T<:Number, TZ<:Number, TP<:Number}
    # Validate input and output dimensions match
    if !(size(z, 1, 2) == size(p, 1, 2) == size(k, 1, 2))
        error("Dimensions for s, p, and k must match")
    end

    TR = promote_type(TZ,TP) # FIXME: Include Complex128
    matrix = Array{SisoZpk{T, TR}}(size(z,1), size(z,2)) # TODO: make this nicer
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
    TR = promote_type(Complex128,T)
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

zpk(sys::StateSpace) = zpk(convert(TransferFunction{SisoRational}, sys))

# We can neither guarantee
function zpk(G::TransferFunction{S}) where {T1, S<:SisoTf{T1}}
    Tnew = typeof(one(T1)/one(T1))
    convert(TransferFunction{SisoZpk{Tnew, complex(Tnew)}}, G)
end

zpk(var::AbstractString) = zpk(tf(var))
zpk(var::AbstractString, Ts::Real) = zpk(tf(var, Ts))
