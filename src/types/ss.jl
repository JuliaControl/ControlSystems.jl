# Convenience constructor for creating StateSpace objects

"""`sys = ss(A, B, C, D, Ts=0) -> sys`


Create a state-space model `sys::StateSpace{T, MT<:AbstractMatrix{T}}`
where `MT` is the type of matrixes `A,B,C,D` and `T` the element type.

This is a continuous-time model if Ts is omitted or set to 0.
Otherwise, this is a discrete-time model with sampling period Ts.
Set Ts=-1 for a discrete-time model with unspecified sampling period.

`sys = ss(D[, Ts, ...])` specifies a static gain matrix D."""
function ss(A::AbstractArray, B::AbstractArray, C::AbstractArray, D::AbstractArray, Ts::Real=0)
    # Check the kwargs for metadata
    nu = size(B, 2)
    ny, nx = size(C, 1), size(C, 2)
    return StateSpace(A, B, C, D, Ts)
end

# Function for accepting scalars
function ss(A::Union{Number,AbstractArray}, B::Union{Number,AbstractArray}, C::Union{Number,AbstractArray}, D::Union{Number,AbstractArray}, Ts::Real=0)
    T = promote_type(eltype(A),eltype(B),eltype(C),eltype(D))
    A = to_matrix(T, A)
    B = to_matrix(T, B)
    C = to_matrix(T, C)
    if D == 0
        D = fill(zero(T), size(C,1), size(B,2))
    else
        D = to_matrix(T, D)
    end
    ss(A, B, C, D, Ts)
end

function ss(A::AbstractArray, B1::AbstractArray, B2::AbstractArray, C1::AbstractArray, C2::AbstractArray, D11::AbstractArray, D12::AbstractArray, D21::AbstractArray, D22::AbstractArray, Ts::Real=0)
    return ExtendedStateSpace(A, B1, B2, C1, C2, D11, D12, D21, D22, Ts)
end

# Function for creation of static gain
function ss(D::AbstractArray{T}, Ts::Real=0) where {T<:Number}
    ny, nu = size(D, 1), size(D, 2)
    A = fill(zero(T), 0, 0)
    B = fill(zero(T), 0, nu)
    C = fill(zero(T), ny, 0)

    return ss(A, B, C, D, Ts)
end
ss(d::Number, Ts::Real=0; kwargs...) = ss([d], Ts)

# ss(sys) converts to StateSpace
ss(sys::LTISystem) = convert(StateSpace, sys)
