"""
    lqr(sys, Q, R)
    lqr(Continuous, A, B, Q, R, args...; kwargs...)
    lqr(Discrete, A, B, Q, R, args...; kwargs...)

Calculate the optimal gain matrix `K` for the state-feedback law `u = -K*x` that
minimizes the cost function:

J = integral(x'Qx + u'Ru, 0, inf) for the continuous-time model `dx = Ax + Bu`.
J = sum(x'Qx + u'Ru, 0, inf) for the discrete-time model `x[k+1] = Ax[k] + Bu[k]`.

Solve the LQR problem for state-space system `sys`. Works for both discrete
and continuous time systems.

The `args...; kwargs...` are sent to the Riccati solver, allowing specification of cross-covariance etc. See `?MatrixEquations.arec / ared` for more help.

To obtain also the solution to the Riccati equation and the eigenvalues of the closed-loop system as well, call `ControlSystemsBase.MatrixEquations.arec / ared` instead (note the different order of the arguments to these functions).

# Examples
Continuous time
```julia
using LinearAlgebra # For identity matrix I
using Plots
A = [0 1; 0 0]
B = [0; 1]
C = [1 0]
sys = ss(A,B,C,0)
Q = I
R = I
L = lqr(sys,Q,R) # lqr(Continuous,A,B,Q,R) can also be used

u(x,t) = -L*x # Form control law,
t=0:0.1:5
x0 = [1,0]
y, t, x, uout = lsim(sys,u,t,x0=x0)
plot(t,x', lab=["Position" "Velocity"], xlabel="Time [s]")
```

Discrete time
```julia
using LinearAlgebra # For identity matrix I
using Plots
Ts = 0.1
A = [1 Ts; 0 1]
B = [0;1]
C = [1 0]
sys = ss(A, B, C, 0, Ts)
Q = I
R = I
L = lqr(Discrete, A,B,Q,R) # lqr(sys,Q,R) can also be used

u(x,t) = -L*x # Form control law,
t=0:Ts:5
x0 = [1,0]
y, t, x, uout = lsim(sys,u,t,x0=x0)
plot(t,x', lab=["Position"  "Velocity"], xlabel="Time [s]")
```
"""
function lqr(::ContinuousType, A, B, Q, R, args...; kwargs...)
    S, _, K = arec(A, B, R, Q, args...; kwargs...)
    return K
end

function lqr(::DiscreteType, A, B, Q, R, args...; kwargs...)
    S, _, K = ared(A, B, R, Q, args...; kwargs...)
    return K
end

@deprecate lqr(A::AbstractMatrix, args...; kwargs...)  lqr(Continuous, A, args...; kwargs...)
@deprecate dlqr(args...; kwargs...)  lqr(Discrete, args...; kwargs...)



"""
    kalman(Continuous, A, C, R1, R2)
    kalman(Discrete, A, C, R1, R2)
    kalman(sys, R1, R2)

Calculate the optimal Kalman gain

The `args...; kwargs...` are sent to the Riccati solver, allowing specification of cross-covariance etc. See `?MatrixEquations.arec/ared` for more help.
"""
kalman(te, A, C, R1,R2, args...; kwargs...) = Matrix(lqr(te, A',C',R1,R2, args...; kwargs...)')

function lqr(sys::AbstractStateSpace, Q, R, args...; kwargs...)
    return lqr(sys.timeevol, sys.A, sys.B, Q, R, args...; kwargs...)
end

function kalman(sys::AbstractStateSpace, R1, R2, args...; kwargs...)
    return Matrix(lqr(sys.timeevol, sys.A', sys.C', R1,R2, args...; kwargs...)')
end

@deprecate kalman(A::AbstractMatrix, args...; kwargs...)  kalman(Continuous, A, args...; kwargs...)
@deprecate dkalman(args...; kwargs...)  kalman(Discrete, args...; kwargs...)

"""
    place(A, B, p, opt=:c)
    place(sys::StateSpace, p, opt=:c)

Calculate the gain matrix `K` such that `A - BK` has eigenvalues `p`.

    place(A, C, p, opt=:o)
    place(sys::StateSpace, p, opt=:o)

Calculate the observer gain matrix `L` such that `A - LC` has eigenvalues `p`.

Uses Ackermann's formula.

**Currently handles only SISO systems**, but a trick is possible to make it work for MIMO systems:
The code below introduces a random projection matrix `P` that projects the inuput space to one dimension, and then shifts the application of `P` from `B` to `K`. 
```julia
nx = 5
nu = 2
A = randn(nx,nx)
B = randn(nx,nu)
P = randn(nu,1)
K = place(A,B*P,zeros(nx))
K2 = P*K
eigvals(A-B*K2)
```

Please note that this function can be numerically sensitive, solving the placement problem in extended precision might be beneficial.
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
