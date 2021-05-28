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
plot(t,x', lab=["Position" "Velocity"], xlabel="Time [s]")
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

function lqr(sys::AbstractStateSpace, Q, R)
    if iscontinuous(sys)
        return lqr(sys.A, sys.B, Q, R)
    else
        return dlqr(sys.A, sys.B, Q, R)
    end
end

function kalman(sys::AbstractStateSpace, R1,R2)
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
plot(t,x', lab=["Position"  "Velocity"], xlabel="Time [s]")
```
"""
function dlqr(A, B, Q, R)
    S = dare(A, B, Q, R)
    K = (B'*S*B + R)\(B'S*A)
    return K
end

function dlqr(sys::AbstractStateSpace, Q, R)
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
    place(A, B, p, opt=:c)
    place(sys::StateSpace, p, opt=:c)

Calculate the gain matrix `K` such that `A - BK` has eigenvalues `p`.

    place(A, C, p, opt=:o)
    place(sys::StateSpace, p, opt=:o)

Calculate the observer gain matrix `L` such that `A - LC` has eigenvalues `p`.

Uses Ackermann's formula.
Currently handles only SISO systems.
"""
function place(A, B, p, opt=:c)
    n = length(p)
    n != size(A,1) && error("Must specify as many poles as states")
    if opt === :c
        n != size(B,1) && error("A and B must have same number of rows")
        if size(B,2) == 1
            acker(A, B, p)
        else
            error("place only implemented for SISO systems")
        end
    elseif opt === :o
        C = B # B is really the "C matrix"
        n != size(C,2) && error("A and C must have same number of columns")
        if size(C,1) == 1
            acker(A', C', p)'
        else
            error("place only implemented for SISO systems")
        end
    else
        error("fourth argument must be :c or :o")
    end
end
function place(sys::AbstractStateSpace, p, opt=:c)
    if opt === :c
        return place(sys.A, sys.B, p, opt)
    elseif opt === :o
        return place(sys.A, sys.C, p, opt)
    else
        error("third argument must be :c or :o")
    end
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
