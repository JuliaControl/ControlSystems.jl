# Convenience constructor for creating StateSpace objects

@doc """`ss(A,B,C,D[, Ts, statenames=..., inputnames=..., outputnames=...]) -> sys`

Create a state-space model.
This is a continuous-time model if Ts is omitted or set to 0.
Otherwise, this is a discrete-time model with sampling period Ts.
Set Ts=-1 for a discrete-time model with unspecified sampling period.

State, input and output names: each can be either a vector of strings (one string per dimension),
or a single string (e.g., "x"). In the latter case, an index is automatically appended to identify
the coordinates for each dimension (e.g. "x1", "x2", ...).

`sys = ss(D[, Ts, ...])` specifies a static gain matrix D.""" ->
function ss(A::Array, B::Array, C::Array, D::Array, Ts::Real=0)
    # Check the kwargs for metadata
    nu = size(B, 2)
    ny, nx = size(C, 1, 2)
    return StateSpace(A, B, C, D, Ts)
end

# Function for accepting scalars
function ss(A::Union{Real,Array}, B::Union{Real,Array}, C::Union{Real,Array}, D::Union{Real,Array}, Ts::Real=0)
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

# Function for creation of static gain
function ss(D::Array{T}, Ts::Real=0) where {T<:Number}
    ny, nu = size(D, 1, 2)
    A = fill(zero(T), 0, 0)
    B = fill(zero(T), 0, nu)
    C = fill(zero(T), ny, 0)

    return ss(A, B, C, D, Ts)
end
ss(d::Real, Ts::Real=0; kwargs...) = ss([d], Ts)

# ss(sys) converts to StateSpace
ss(sys::LTISystem) = convert(StateSpace, sys)
