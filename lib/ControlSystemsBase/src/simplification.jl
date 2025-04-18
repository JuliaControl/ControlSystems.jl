"""
    sminreal(sys)

Compute the structurally minimal realization of the state-space system `sys`. A
structurally minimal realization is one where only states that can be
determined to be uncontrollable and unobservable based on the location of 0s in
`sys` are removed.

Systems with numerical noise in the coefficients, e.g., noise on the order of `eps` require truncation to zero to be
affected by structural simplification, e.g.,
```julia
trunc_zero!(A) = A[abs.(A) .< 10eps(maximum(abs, A))] .= 0
trunc_zero!(sys.A); trunc_zero!(sys.B); trunc_zero!(sys.C)
sminreal(sys)
```

In contrast to [`minreal`](@ref), which performs pole-zero cancellation using linear-algebra operations, has an 𝑂(nₓ^3) complexity and is subject to numerical tolerances, `sminreal` is computationally very cheap and numerically exact (operates on integers). However, the ability of `sminreal` to reduce the order of the model is much less powerful.

See also [`minreal`](@ref).
"""
function sminreal(sys::AbstractStateSpace)
    A, B, C, inds = struct_ctrb_obsv(sys)
    return basetype(sys)(A, B, C, sys.D, sys.timeevol)
end

# Determine the structurally controllable and observable realization for the system
struct_ctrb_obsv(sys::AbstractStateSpace) = struct_ctrb_obsv(sys.A, sys.B, sys.C)

function struct_ctrb_obsv(A::AbstractVecOrMat, B::AbstractVecOrMat, C::AbstractVecOrMat)
    costates = struct_ctrb_states(A, B) .& struct_ctrb_states(A', C')
    if !all(costates)
        inds = findall(costates)
        return A[inds, inds], B[inds, :], C[:, inds], inds
    else
        return A, B, C, [1:size(A, 1);]
    end
end

"""
    struct_ctrb_states(A::AbstractVecOrMat, B::AbstractVecOrMat)

Compute a bit-vector, expressing whether a state of the pair (A, B) is
structurally controllable based on the location of zeros in the matrices. 
"""
function struct_ctrb_states(A::AbstractVecOrMat, B::AbstractVecOrMat)
    size(A,1) > typemax(UInt16) && error("Maximum size of A exceeded. If you encounter this error, please open an issue. This limit is not fundamental and exists for performance reasons only.")
    # UInt16 can only store up to 65535, so if A is completely dense and of size larger than 65535, the computations below might overflow. This is exceedingly unlikely though.
    bitA = UInt16.(.!iszero.(A)) # Convert to Int because multiplying with a bit matrix is slow
    x = vec(any(B .!= 0, dims=2)) # index vector indicating states that have been affected by input
    if A isa SMatrix
        xi = Vector(bitA * x) # To aboid setindex! error below
    else
        xi = bitA * x
    end
    xi2 = similar(xi)
    @. xi = (xi != false) | !iszero(x)
    for i = 2:size(A, 1) # apply A nx times, similar to controllability matrix
        mul!(xi2, bitA, xi)
        @. xi = (xi2 != false) | !iszero(xi)
    end
    xi .!= false # Convert back to BitVector
end