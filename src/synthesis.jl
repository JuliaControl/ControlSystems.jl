@doc """`lqr(A, B, Q, R)`

Calculate the optimal gain matrix `K` for the state-feedback law `u = K*x` that
minimizes the cost function:

J = integral(x'Qx + u'Ru, 0, inf).

For the continuous time model `dx = Ax + Bu`.

`lqr(sys, Q, R)`

Solve the LQR problem for state-space system `sys`. Works for both discrete
and continuous time systems.""" ->
function lqr{T<:BlasFloat}(A::StridedMatrix{T}, B::StridedMatrix{T},
            Q::StridedMatrix{T}, R::StridedMatrix{T})
    S = care(A, B, Q, R)
    K = R\B'*S
    return K
end
lqr{T<:Integer}(A::StridedMatrix{T}, B::StridedMatrix{T}, Q::StridedMatrix{T},
        R::StridedMatrix{T}) = lqr(float(A), float(B), float(Q), float(R))

function lqr(sys::StateSpace, Q, R)
    if iscontinuous(sys)
        return lqr(sys.A, sys.B, Q, R)
    else
        return dlqr(sys.A, sys.B, Q, R)
    end
end

@doc """`dlqr(A, B, Q, R)`, `dlqr(sys, Q, R)`

Calculate the optimal gain matrix `K` for the state-feedback law `u[k] = K*x[k]` that
minimizes the cost function:

J = sum(x'Qx + u'Ru, 0, inf).

For the discrte time model `x[k+1] = Ax[k] + Bu[k]`.""" ->
function dlqr{T<:BlasFloat}(A::StridedMatrix{T}, B::StridedMatrix{T},
            Q::StridedMatrix{T}, R::StridedMatrix{T})
    S = dare(A, B, Q, R)
    K = (B'*S*B + R)\(B'S*A)
    return K
end
dlqr{T<:Integer}(A::StridedMatrix{T}, B::StridedMatrix{T}, Q::StridedMatrix{T},
        R::StridedMatrix{T}) = dlqr(float(A), float(B), float(Q), float(R))
