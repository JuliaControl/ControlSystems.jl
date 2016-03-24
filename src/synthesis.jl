@doc """`lqr(A, B, Q, R)`

Calculate the optimal gain matrix `K` for the state-feedback law `u = K*x` that
minimizes the cost function:

J = integral(x'Qx + u'Ru, 0, inf).

For the continuous time model `dx = Ax + Bu`.

`lqr(sys, Q, R)`

Solve the LQR problem for state-space system `sys`. Works for both discrete
and continuous time systems.""" ->
function lqr(A, B, Q, R)
    S = care(A, B, Q, R)
    K = R\B'*S
    return K
end

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
function dlqr(A, B, Q, R)
    S = dare(A, B, Q, R)
    K = (B'*S*B + R)\(B'S*A)
    return K
end

@doc """`place(A, B, p)`, `place(sys::StateSpace, p)`

Calculate gain matrix `K` such that
the poles of `(A-BK)` in are in `p`""" ->
function place(A, B, p)
    n = length(p)
    n != size(A,1) && error("Must define as many poles as states")
    n != size(B,1) && error("A and B must have same number of rows")
    if size(B,2) == 1
        acker(A,B,p)
    else
        error("place only implemented for SISO systems")
    end
end

function place(sys::StateSpace, p)
    return place(sys.A, sys.B, p)
end

#Implements Ackermann's formula for placing poles of (A-BK) in p
function acker(A,B,P)
  n = length(P)
  #Calculate characteristic polynomial
  poly = reduce(*,Poly([1]),[Poly([1, -p]) for p in P])
  q = zero(Array{promote_type(eltype(A),Float64),2}(n,n))
  for i = n:-1:0
      q += A^(n-i)*poly[i+1]
  end
  S = Array{promote_type(eltype(A),eltype(B),Float64),2}(n,n)
  for i = 0:(n-1)
      S[:,i+1] = A^i*B
  end
  return [zeros(1,n-1) 1]*(S\q)
end


"""
`feedback(L)` Return L/(1+L)
`feedback(P,C)` Return PC/(1+PC)
"""
feedback(L::TransferFunction) = L/(1+L)
feedback(P::TransferFunction, C::TransferFunction) = feedback(P*C)


"""
`feedback2dof(P,R,S,T)` Return `BT/(AR+ST)` where B and A are the numerator and denomenator polynomials of `P` respectively
`feedback2dof(B,A,R,S,T)` Return `BT/(AR+ST)`
"""
feedback2dof(P::TransferFunction,R,S,T) = tf(conv(tfnum(P),T),zpconv(tfden(P),R,tfnum(P),S))
feedback2dof(B,A,R,S,T) = tf(conv(B,T),zpconv(A,R,B,S))
