"""
    lqr(A, B, Q, R)

Calculate the optimal gain matrix `K` for the state-feedback law `u = -K*x` that
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

u(x,t) = -L*x # Form control law,
t=0:0.1:5
x0 = [1,0]
y, t, x, uout = lsim(sys,u,t,x0=x0)
plot(t,x, lab=["Position" "Velocity"], xlabel="Time [s]")
```
"""
function lqr(A, B, Q, R)
    S = care(A, B, Q, R)
    K = R\B'*S
    return K
end

"""
    kalman(A, C, R1, R2)
    kalman(sys, R1, R2)

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


"""
    dlqr(A, B, Q, R)
    dlqr(sys, Q, R)

Calculate the optimal gain matrix `K` for the state-feedback law `u[k] = -K*x[k]` that
minimizes the cost function:

J = sum(x'Qx + u'Ru, 0, inf).

For the discrte time model `x[k+1] = Ax[k] + Bu[k]`.

See also `lqg`

Usage example:
```julia
using LinearAlgebra # For identity matrix I
Ts = 0.1
A = [1 Ts; 0 1]
B = [0;1]
C = [1 0]
sys = ss(A,B,C,0, Ts)
Q = I
R = I
L = dlqr(A,B,Q,R) # lqr(sys,Q,R) can also be used

u(x,t) = -L*x # Form control law,
t=0:Ts:5
x0 = [1,0]
y, t, x, uout = lsim(sys,u,t,x0=x0)
plot(t,x, lab=["Position"  "Velocity"], xlabel="Time [s]")
```
"""
function dlqr(A, B, Q, R)
    S = dare(A, B, Q, R)
    K = (B'*S*B + R)\(B'S*A)
    return K
end

function dlqr(sys::StateSpace, Q, R)
    !isdiscrete(sys) && throw(ArgumentError("Input argument sys must be discrete-time system"))
    return dlqr(sys.A, sys.B, Q, R)
end

"""
    dkalman(A, C, R1, R2)
    dkalman(sys, R1, R2)

Calculate the optimal Kalman gain for discrete time systems

"""
dkalman(A, C, R1,R2) = Matrix(dlqr(A',C',R1,R2)')

"""
    place(A, B, p)
    place(sys::StateSpace, p)

Calculate gain matrix `K` such that
the poles of `(A-BK)` in are in `p`.

Uses Ackermann's formula.
For observer pole placement, see `luenberger`.
"""
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

"""
    luenberger(A, C, p)
    luenberger(sys::StateSpace, p)

Calculate gain matrix `L` such that the poles of `(A - LC)` are in `p`.
Uses sytem's dual form (Controllability-Observability duality) applied to Ackermann's formula.
That is, `(A - BK)` is indentic to `(A' - C'L') == (A - LC)`.
"""
function luenberger(A, C, p)
    place(A', C', p)'
end

function luenberger(sys::StateSpace, p)
    return luenberger(sys.A, sys.C, p)
end

#Implements Ackermann's formula for placing poles of (A-BK) in p
function acker(A,B,P)
    n = length(P)
    #Calculate characteristic polynomial
    poly = mapreduce(p -> Polynomial([1, -p]), *, P, init=Polynomial(one(eltype(P))))
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
