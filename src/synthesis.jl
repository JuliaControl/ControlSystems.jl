"""`lqr(A, B, Q, R)`

Calculate the optimal gain matrix `K` for the state-feedback law `u = K*x` that
minimizes the cost function:

J = integral(x'Qx + u'Ru, 0, inf).

For the continuous time model `dx = Ax + Bu`.

`lqr(sys, Q, R)`

Solve the LQR problem for state-space system `sys`. Works for both discrete
and continuous time systems.

See also `LQG`

Usage example:
```julia
using LinearAlgebra # For identity matrix I
A = [0 1; 0 0]
B = [0;1]
C = [1 0]
sys = ss(A,B,C,0)
Q = I
R = I
L = lqr(sys,Q,R)

u(t,x) = -L*x # Form control law,
t=0:0.1:5
x0 = [1,0]
y, t, x, uout = lsim(sys,u,t,x0)
plot(t,x, lab=["Position", "Velocity"]', xlabel="Time [s]")
```
"""
function lqr(A, B, Q, R)
    S = care(A, B, Q, R)
    K = R\B'*S
    return K
end

"""`kalman(A, C, R1, R2)` kalman(sys, R1, R2)`

Calculate the optimal Kalman gain

See also `LQG`

"""
kalman(A, C, R1,R2) = Matrix(lqr(A',C',R1,R2)')

function lqr(sys::StateSpace, Q, R)
    if iscontinuous(sys)
        return lqr(sys.A, sys.B, Q, R)
    else
        return dlqr(sys.A, sys.B, Q, R)
    end
end

function kalman(sys::StateSpace, R1,R2)
    if iscontinuous(sys)
        return Matrix(lqr(sys.A', sys.C', R1,R2)')
    else
        return Matrix(dlqr(sys.A', sys.C', R1,R2)')
    end
end


"""`dlqr(A, B, Q, R)`, `dlqr(sys, Q, R)`

Calculate the optimal gain matrix `K` for the state-feedback law `u[k] = K*x[k]` that
minimizes the cost function:

J = sum(x'Qx + u'Ru, 0, inf).

For the discrte time model `x[k+1] = Ax[k] + Bu[k]`.

See also `lqg`

Usage example:
```julia
using LinearAlgebra # For identity matrix I
h = 0.1
A = [1 h; 0 1]
B = [0;1]
C = [1 0]
sys = ss(A,B,C,0, h)
Q = I
R = I
L = dlqr(A,B,Q,R) # lqr(sys,Q,R) can also be used

u(t,x) = -L*x # Form control law,
t=0:h:5
x0 = [1,0]
y, t, x, uout = lsim(sys,u,t,x0)
plot(t,x, lab=["Position", "Velocity"]', xlabel="Time [s]")
```
"""
function dlqr(A, B, Q, R)
    S = dare(A, B, Q, R)
    K = (B'*S*B + R)\(B'S*A)
    return K
end

function dlqr(A, B::Vector, Q, R)
    dlqr(A, reshape(B, length(B), 1), Q, R)
end

"""`dkalman(A, C, R1, R2)` kalman(sys, R1, R2)`

Calculate the optimal Kalman gain for discrete time systems

"""
dkalman(A, C, R1,R2) = Matrix(dlqr(A',C',R1,R2)')

"""`place(A, B, p)`, `place(sys::StateSpace, p)`

Calculate gain matrix `K` such that
the poles of `(A-BK)` in are in `p`.

Uses Ackermann's formula."""
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
    poly = mapreduce(p -> Poly([1, -p]), *, P, init=Poly(one(eltype(P))))
    q = zero(Array{promote_type(eltype(A),Float64),2}(undef, n,n))
    for i = n:-1:0
        q += A^(n-i)*poly[i]
    end
    S = Array{promote_type(eltype(A),eltype(B),Float64),2}(undef, n,n)
    for i = 0:(n-1)
        S[:,i+1] = A^i*B
    end
    return [zeros(1,n-1) 1]*(S\q)
end


"""
`feedback(L)` Returns L/(1+L)
`feedback(P1,P2)` Returns P1/(1+P1*P2)
"""
feedback(L::TransferFunction) = L/(1+L)
feedback(P1::TransferFunction, P2::TransferFunction) = P1/(1+P1*P2)

#Efficient implementations
function feedback(L::TransferFunction{T}) where T<:SisoRational
    if size(L) != (1,1)
        error("MIMO TransferFunction feedback isn't implemented, use L/(1+L)")
    end
    P = numpoly(L)
    Q = denpoly(L)
    tf(P, P+Q, L.Ts)
end

function feedback(L::TransferFunction{T}) where {T<:SisoZpk}
    if size(L) != (1,1)
        error("MIMO TransferFunction feedback isn't implemented, use L/(1+L)")
    end
    #Extract polynomials and create P/(P+Q)
    k = L.matrix[1].k
    denpol = numpoly(L)[1]+denpoly(L)[1]
    kden = denpol[end] # Get coeff of s^n
    # Create siso system
    sisozpk = T(L.matrix[1].z, roots(denpol), k/kden)
    return TransferFunction{T}(fill(sisozpk,1,1), L.Ts)
end

"""
`feedback(sys)`

`feedback(sys1,sys2)`

Forms the negative feedback interconnection
```julia
>-+ sys1 +-->
  |      |
 (-)sys2 +
```
If no second system is given, negative identity feedback is assumed
"""
function feedback(sys::StateSpace)
    sys.ny != sys.nu && error("Use feedback(sys1::StateSpace,sys2::StateSpace) if sys.ny != sys.nu")
    feedback(sys,ss(Matrix{numeric_type(sys)}(I,sys.ny,sys.ny)))
end

function feedback(sys1::StateSpace,sys2::StateSpace)
    sum(abs.(sys1.D)) != 0 && sum(abs.(sys2.D)) != 0 && error("There can not be a direct term (D) in both sys1 and sys2")
    A = [sys1.A+sys1.B*(-sys2.D)*sys1.C sys1.B*(-sys2.C); sys2.B*sys1.C  sys2.A+sys2.B*sys1.D*(-sys2.C)]
    B = [sys1.B; sys2.B*sys1.D]
    C = [sys1.C  sys1.D*(-sys2.C)]
    ss(A,B,C,sys1.D)
end


"""
`feedback2dof(P,R,S,T)` Return `BT/(AR+ST)` where B and A are the numerator and denomenator polynomials of `P` respectively
`feedback2dof(B,A,R,S,T)` Return `BT/(AR+ST)`
"""
function feedback2dof(P::TransferFunction,R,S,T)
    !issiso(P) && error("Feedback not implemented for MIMO systems")
    tf(conv(poly2vec(numpoly(P)[1]),T),zpconv(poly2vec(denpoly(P)[1]),R,poly2vec(numpoly(P)[1]),S))
end

feedback2dof(B,A,R,S,T) = tf(conv(B,T),zpconv(A,R,B,S))
