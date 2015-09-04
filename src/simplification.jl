@doc """`sminreal(sys)`

Compute the structurally minimal realization of the state-space system `sys`. A
structurally minimal realization is one where only states that can be
determined to be uncontrollable and unobservable based on the location of 0s in
`sys` are removed.""" ->
function sminreal(sys::StateSpace)
    A, B, C, inds = struct_ctrb_obsv(sys)
    return StateSpace(A, B, C, sys.D, sys.Ts, sys.statenames[inds],
            sys.inputnames, sys.outputnames)
end

# Determine the structurally controllable and observable realization for the system
struct_ctrb_obsv(sys::StateSpace) = struct_ctrb_obsv(sys.A, sys.B, sys.C)
function struct_ctrb_obsv(A::VecOrMat, B::VecOrMat, C::VecOrMat)
    costates = struct_ctrb_states(A, B) & struct_ctrb_states(A', C')
    if !all(costates)
        inds = find(costates)
        return A[inds, inds], B[inds, :], C[:, inds], inds
    else
        return A, B, C, [1:size(A, 1);]
    end
end

# Compute a bit-vector, expressing whether a state of the pair (A, B) is
# structurally controllable.
function struct_ctrb_states(A::VecOrMat, B::VecOrMat)
    bitA = A .!= 0
    d_cvec = cvec = any(B .!= 0, 2)
    while any(d_cvec .!= 0)
        Adcvec = any(bitA[:, find(d_cvec)], 2)
        cvec = cvec | Adcvec
        d_cvec = Adcvec & ~cvec
    end
    return cvec
end
