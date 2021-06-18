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
"""
function sminreal(sys::StateSpace)
    A, B, C, inds = struct_ctrb_obsv(sys)
    return StateSpace(A, B, C, sys.D, sys.timeevol)
end

# Determine the structurally controllable and observable realization for the system
struct_ctrb_obsv(sys::StateSpace) = struct_ctrb_obsv(sys.A, sys.B, sys.C)

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
    bitA = A .!= 0
    x = vec(any(B .!= 0, dims=2)) # index vector indicating states that have been affected by input
    for i = 1:size(A, 1) # apply A nx times, similar to controllability matrix
        x = x .| (bitA * x) .!= 0
    end
    x
end
