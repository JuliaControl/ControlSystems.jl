function tf(num::AbstractVecOrMat{<:AbstractVector{T1}}, den::AbstractVecOrMat{<:AbstractVector{T2}}, Ts::Real=0.0) where {T1<:Number, T2<:Number}
    # Validate input and output dimensions match
    ny, nu = size(num, 1, 2)
    if (ny, nu) != size(den, 1, 2)
        error("num and den dimensions must match")
    end

    T = promote_type(T1,T2)
    matrix = Matrix{SisoRational{T}}(ny, nu)
    for o=1:ny
        for i=1:nu
            matrix[o, i] = SisoRational{T}(Vector{T}(num[o, i]), Vector{T}(den[o, i]))
        end
    end
    return TransferFunction{SisoRational{T}}(matrix, Ts)
end
function tf(num::AbstractVector{T1}, den::AbstractVector{T2}, Ts::Real=0.0) where {T1<:Number, T2<:Number}
    T = promote_type(T1, T2)
    return TransferFunction{SisoRational{T}}(fill(SisoRational{T}(num, den), 1, 1), Ts)
end
tf(num::Real, den::Vector, Ts::Real=0.0) = tf([num], den, Ts)

# Cases for just static gain
function tf(D::AbstractArray{T}, Ts::Real=0.0) where {T<:Number}
    ny, nu = size(D, 1, 2)

    matrix = Matrix{SisoRational{T}}(ny, nu)
    for i in eachindex(D)
        matrix[i] = SisoRational{T}([D[i]], [one(T)])
    end
    return TransferFunction{SisoRational{T}}(matrix, Float64(Ts))
end
tf(n::Real, Ts::Real=0; kwargs...) = tf([n], Ts; kwargs...)

# Function for creation of 's' or 'z' var
function tf(var::AbstractString)
    var != "s" && error("var must be 's' for continuous time tf.")
    return tf([1, 0], [1])
end
function tf(var::AbstractString, Ts::Real)
    var != "z" && error("var must be 'z' for discrete time tf.")
    Ts == 0 && error("Ts must not be 0 for discrete time tf.")
    return tf([1, 0], [1], Ts)
end
