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

u(x,t) = -L*x # Form control law,
t=0:h:5
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

"""`dkalman(A, C, R1, R2)` kalman(sys, R1, R2)`

Calculate the optimal Kalman gain for discrete time systems

"""
dkalman(A, C, R1,R2) = Matrix(dlqr(A',C',R1,R2)')

"""`place(A, B, p)`, `place(sys::StateSpace, p)`

Calculate gain matrix `K` such that
the poles of `(A-BK)` in are in `p`.

Uses Ackermann's formula."""
function place(A, B, p; kwargs...)
    n = length(p)
    n != size(A,1) && error("Must define as many poles as states")
    n != size(B,1) && error("A and B must have same number of rows")
    if size(B,2) == 1
        acker(A,B,p)
    else 
        placemimo(A, B, p; kwargs...)
    end
end

function place(sys::StateSpace, p; kwargs...)
    return place(sys.A, sys.B, p; kwargs...)
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

# Implemented according to method 0 in https://perso.uclouvain.be/paul.vandooren/publications/KautskyNV85.pdf
# The method is not guaranteed to converge (compared to method 1) but is much simpler and faster.
# A and B should be real and controllable, p is complex and closed under conjugation
function placemimo(A, B, p; max_iter=10000, tol=1e-7)
    n, m = size(B)

    U, Z = qr(B) # Assume B is full rank. Should this be assertion?
    
    # Check for more errors?

    # Check for controllability (just temporary)
    #@assert isposdef(gram(ss(A, B, zeros(1, size(A, 2)), zeros(1, size(B, 2))), :c)) "(A, B) is not controllable"
    @assert rank(ctrb(A, B)) == size(A, 1) "(A, B) is not controllable"

    if n <= m # If n <= m we use X=I and just need to calculate the last step
        return -Z \ U' * (Diagonal(p) - A) 
    else # We find an initial X and iteratively improve it
        U₀, U₁ = U[:, 1:m], U[:, m+1:end]
        # TODO is this type correct? Probably not
        T = promote_type(eltype(U), eltype(p))
        S = Vector{Matrix{T}}(undef, length(p))
        for i in eachindex(p)
            Q, R = qr((A' - p[i]*I) * U₁) # Here we don't assume full rank and need to check where to cut Q
            S[i] = Q[:, rank(R)+1:end] # Extract vectors describing nullspace of U1'*(A-pI), can we do this more efficient? rank(R) = m if (A, B) controllable
        end

        X = Matrix{T}(undef, size(A))
        for j in 1:n
            proj = adjoint(S[j]) * randn(T, n, 1) # Should random here be complex if complex poles?
            X[:, j] = S[j] * proj / norm(proj)
        end
        # Think the init should work like this with high prob, otherwise could do xi = Si*yi and find yi with some linear optimization?
        # Should this work as well with identity matrix projected on the nullspaces?
        # Seems like this will be a problem when #same poles > m since m sets the dimension of the nullspaces Sj and same pole => same Sj
        @assert rank(X) == n "X=$X is incorrectly initialized"
        
        ν₂ = cond(X)
        for i in 1:max_iter
            for j in 1:n
                Q, = qr(X[:, 1:end .!= j])
                γ = Q[:, end]                       # γ is the most orthogonal vector to [x1, x2, ..., xj-1, xj+1, ..., xn]
                proj = adjoint(S[j]) * γ                    # project γ onto the nullspace for this element, and set xj to that
                X[:, j] = S[j] * proj / norm(proj)
            end

            nν₂ = cond(X)
            #@show ν₂, nν₂
            if ν₂ - nν₂ < tol  
                M = transpose(lu(transpose(X)) \ transpose(X * Diagonal(p)))
                return -Z \ U₀' * (M - A) 
            end
            #ν₂ = nν₂
            ν₂ = (nν₂ + ν₂) / 2
        end
        error("Condition number did not converge, try to increase `max_iter` or reduce `tol`.")
    end
end