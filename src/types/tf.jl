""" `sys = tf(num, den[, Ts]), sys = tf(gain[, Ts])`

Create as a fraction of polynomials:

`sys::TransferFunction{SisoRational{T,TR}} = numerator/denominator`
where T is the type of the coefficients in the polynomial.

`num`: the coefficients of the numerator polynomial. Either scalar or vector to create SISO systems
or an array of vectors to create MIMO system.

`den`: the coefficients of the denominator polynomial. Either vector to create SISO systems
or an array of vectors to create MIMO system.

`Ts`: Sample time if discrete time system.

Other uses:
`tf(sys)`: Convert `sys` to `tf` form.
`tf("s")`, `tf("z")`: Create the continuous transferfunction `s`.

See also: `zpk`, `ss`
"""
function tf(num::AbstractVecOrMat{<:AbstractVector{T1}}, den::AbstractVecOrMat{<:AbstractVector{T2}}, Ts::TE) where {TE<:TimeEvolution, T1<:Number, T2<:Number}
    # Validate input and output dimensions match
    ny, nu = size(num, 1), size(num, 2)
    if (ny, nu) != (size(den, 1), size(den, 2))
        error("num and den dimensions must match")
    end

    T = promote_type(T1,T2)
    matrix = Matrix{SisoRational{T}}(undef, ny, nu)
    for o=1:ny
        for i=1:nu
            matrix[o, i] = SisoRational{T}(Vector{T}(num[o, i]), Vector{T}(den[o, i]))
        end
    end
    return TransferFunction{TE,SisoRational{T}}(matrix, Ts)
end
tf(num::AbstractVecOrMat{<:AbstractVector{T1}}, den::AbstractVecOrMat{<:AbstractVector{T2}}, Ts::Number) where {T1,T2} =
    tf(num, den, Discrete(Ts))
tf(num::AbstractVecOrMat{<:AbstractVector{T1}}, den::AbstractVecOrMat{<:AbstractVector{T2}}) where {T1,T2} =
    tf(num, den, Continuous())

function tf(num::AbstractVector{T1}, den::AbstractVector{T2}, timeevol::TE) where {TE<:TimeEvolution,T1<:Number, T2<:Number}
    T = promote_type(T1, T2)
    return TransferFunction{TE,SisoRational{T}}(fill(SisoRational{T}(num, den), 1, 1), timeevol)
end
tf(num::AbstractVector{T1}, den::AbstractVector{T2}, Ts::Number) where {T1<:Number, T2<:Number} =
    tf(num, den, Discrete(Ts))
tf(num::AbstractVector{T1}, den::AbstractVector{T2}) where {T1<:Number, T2<:Number} =
    tf(num, den, Continuous())

tf(num::Number, den::Vector, args...) = tf([num], den, args...)

# Cases for just static gain
function tf(D::AbstractArray{T}, timeevol::TE) where {TE<:TimeEvolution, T<:Number}
    ny, nu = size(D, 1), size(D, 2)

    matrix = Matrix{SisoRational{T}}(undef, ny, nu)
    for i in eachindex(D)
        matrix[i] = SisoRational{T}([D[i]], [one(T)])
    end
    return TransferFunction{TE,SisoRational{T}}(matrix, timeevol)
end
tf(D::AbstractArray{T}, Ts::Number) where T = tf(D, Discrete(Ts))
tf(D::AbstractArray{T}) where T = tf(D, Continuous())

tf(n::Number, args...; kwargs...) = tf([n], args...; kwargs...)

tf(sys::AbstractStateSpace) = convert(TransferFunction, sys) # NOTE: Would perhaps like to write TransferFunction{SisoRational}, but couldn't get this to work...

function tf(G::TransferFunction{TE,<:SisoTf{T}}) where {TE<:TimeEvolution,T<:Number}
    convert(TransferFunction{TE,SisoRational{T}}, G)
end

# Function for creation of 's' or 'z' var
function tf(var::AbstractString)
    var != "s" && error("var must be 's' for continuous time tf.")
    return tf([1, 0], [1], Continuous())
end
function tf(var::AbstractString, Ts::Real)
    var != "z" && error("var must be 'z' for discrete time tf.")
    Ts == 0 && error("Ts must not be 0 for discrete time tf.")
    return tf([1, 0], [1], Discrete(Ts))
end

## Constructors for polynomial inputs
function tf(num::AbstractArray{PT}, den::AbstractArray{PT},  timeevol::TE) where {TE<:TimeEvolution,T<:Number, PT <: Polynomials.Polynomial{T}}
    ny, nu = size(num, 1), size(num, 2)
    if (ny, nu) != (size(den, 1), size(den, 2))
        error("num and den dimensions must match")
    end

    matrix = Matrix{SisoRational{T}}(undef, ny, nu)
    for o=1:ny
        for i=1:nu
            matrix[o, i] = SisoRational{T}(num[o, i], den[o, i])
        end
    end
    return TransferFunction{TE,SisoRational{T}}(matrix, timeevol)
end
tf(num::AbstractArray{PT}, den::AbstractArray{PT}, Ts::Number) where {T<:Number, PT <: Polynomials.Polynomial{T}} =
    tf(num, den, Discrete(Ts))
tf(num::AbstractArray{PT}, den::AbstractArray{PT}) where {T<:Number, PT <: Polynomials.Polynomial{T}} =
    tf(num, den, Continuous())

function tf(num::PT, den::PT, timeevol::TE) where {TE<:TimeEvolution, T<:Number, PT <: Polynomials.Polynomial{T}}
    tf(fill(num,1,1), fill(den,1,1), timeevol)
end
tf(num::PT, den::PT, Ts::Number) where {T<:Number, PT <: Polynomials.Polynomial{T}} =
    tf(num, den, Discrete(Ts))
tf(num::PT, den::PT) where {T<:Number, PT <: Polynomials.Polynomial{T}} =
    tf(num, den, Continuous())
