"""
    lqr(sys,              Q, R;          extra = Val(false))
    lqr(Continuous, A, B, Q, R, args...; extra = Val(false), kwargs...)
    lqr(Discrete,   A, B, Q, R, args...; extra = Val(false), kwargs...)

Calculate the optimal gain matrix `K` for the state-feedback law `u = -K*x` that
minimizes the cost function:

J = integral(x'Qx + u'Ru, 0, inf) for the continuous-time model `dx = Ax + Bu`.
J = sum(x'Qx + u'Ru, 0, inf) for the discrete-time model `x[k+1] = Ax[k] + Bu[k]`.

Solve the LQR problem for state-space system `sys`. Works for both discrete
and continuous time systems.

The `args...; kwargs...` are sent to the Riccati solver, allowing specification of cross-covariance etc. See `?MatrixEquations.arec / ared` for more help.

To obtain also the solution to the Riccati equation and the eigenvalues of the closed-loop system as well, call `ControlSystemsBase.MatrixEquations.arec / ared` instead (note the different order of the arguments to these functions).

To obtain a discrete-time approximation to a continuous-time LQR problem, the function [`c2d`](@ref) can be used to obtain corresponding discrete-time cost matrices.

If `extra = Val(true)`, the function returns `K, P, p`, where `P` is the solution to the Riccati equation (Lyapunov function `x'P*x`), and `p` are the eigenvalues of `A-BK`.

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

# FAQ
This function requires
- `Q` must be positive semi-definite
- `R` must be positive definite
- The pair `(Q,A)` must not have any unobservable modes on the imaginary axis (cont) / unit circle (disc), e.g., there must not be any integrating modes that are not penalized by `Q`. if this condition does not hold, you may get the error "The Hamiltonian matrix is not dichotomic".
"""
function lqr(::ContinuousType, A, B, Q, R, args...; extra::Val{E} = Val(false), kwargs...) where E
    S, p, K, args... = arec(A, B, R, Q, args...; kwargs...)
    E ? (K, S, p, args...) : K
end

function lqr(::DiscreteType, A, B, Q, R, args...; extra::Val{E} = Val(false), kwargs...) where E
    S, p, K, args... = ared(A, B, R, Q, args...; kwargs...)
    E ? (K, S, p, args...) : K
end

@deprecate lqr(A::AbstractMatrix, args...; kwargs...)  lqr(Continuous, A, args...; kwargs...)
@deprecate dlqr(args...; kwargs...)  lqr(Discrete, args...; kwargs...)



"""
    kalman(Continuous, A, C, R1, R2; extra=Val(false))
    kalman(Discrete,   A, C, R1, R2; extra=Val(false), direct = false)
    kalman(sys,              R1, R2; extra=Val(false), direct = false)

Calculate the optimal asymptotic Kalman gain for the linear-Gaussian model
```math
\\begin{aligned}
dx &= Ax + Bu + w \\\\
y &= Cx + v
\\end{aligned}
```
where `w` is the dynamics noise with covariance `R1` and `v` is the measurement noise with covariance `R2`.

In discrete time, the returned Kalman gain `K` is designed to be used with (`direct = false`)
```math
x(t+1|t) = (A-KC)x(t|t-1) + Bu(t) + Ky(t)
```
and (`direct = true`)
```math
x(t+1|t+1) = (A-KC)x(t|t) + Bu(t) + Ky(t+1)
```

If `direct = true`, the observer gain is computed as ``K_d = R_∞C^T (R_2 + CR_∞ C^T)^{-1}`` instead of ``K = A K_d``. This option is intended to be used together with the option `direct = true` to [`observer_controller`](@ref). Ref: "Computer-Controlled Systems" pp 140. `direct = false` is sometimes referred to as a "delayed" estimator, while `direct = true` is a "current" estimator.

To obtain a discrete-time approximation to a continuous-time LQG problem, the function [`c2d`](@ref) can be used to obtain corresponding discrete-time covariance matrices.

To obtain an LTISystem that represents the Kalman filter, pass the obtained Kalman feedback gain into [`observer_filter`](@ref). To obtain an LQG controller, pass the obtained Kalman feedback gain as well as a state-feedback gain computed using [`lqr`](@ref) into [`observer_controller`](@ref).

If `extra = Val(true)`, the function returns `K, R∞, p`, where `R∞` is the solution to the Riccati equation (the stationary prediction-error covariance matrix, the covariance of ``x̃(t|t-1)``), and `p` are the eigenvalues of `A-KC`. In this case, the _filtering error_ covariance is given by ``(I-KC) R∞``, i.e., the covariance of ``x̃(t|t)``.

The `args...; kwargs...` are sent to the Riccati solver, allowing specification of cross-covariance etc. See `?ControlSystemsBase.MatrixEquations.arec/ared` for more help.

# FAQ
This function requires
- `R1` must be positive semi-definite
- `R2` must be positive definite
- The pair `(A,R1)` must not have any uncontrollable modes on the imaginary axis (cont) / unit circle (disc), e.g., there must not be any integrating modes that are not affected through `R1`. if this condition does not hold, you may get the error "The Hamiltonian matrix is not dichotomic".

# Extended help
`MatrixEquations.ared` solves the Riccati equation corresponding to the direct form, but returns the ``K`` matrix for the indirect form. The solution to the Riccati equation is the stationary prediction-error covariance matrix ``R_∞``, and the filtering error covariance is given by ``(I-A^{-1}KC) R_∞``. If using the direct form, the filter-error covariance is given by ``(I-K_d C) R_∞`` 
"""
function kalman(te, A, C, R1,R2, args...; direct = false, extra::Val{E} = Val(false), kwargs...) where E
    Kt, R∞, p, args... = lqr(te, A',C',R1,R2, args...; extra=Val(true), kwargs...)
    if direct
        te isa ContinuousType && error("direct = true only applies to discrete-time systems")
        K = (R∞*C')/(R2 + C*R∞*C') # ared returns K for the indirect form, which is A*Kdirect
        return E ? (K, R∞, p, args...) : K
    end
    E ? (Matrix(Kt'), R∞, p, args...) : Matrix(Kt')
end

function lqr(sys::AbstractStateSpace, Q, R, args...; kwargs...)
    return lqr(sys.timeevol, sys.A, sys.B, Q, R, args...; kwargs...)
end

function kalman(sys::AbstractStateSpace, R1, R2, args...; kwargs...)
    return kalman(sys.timeevol, sys.A, sys.C, R1,R2, args...; kwargs...)
end

@deprecate kalman(A::AbstractMatrix, args...; kwargs...)  kalman(Continuous, A, args...; kwargs...)
@deprecate dkalman(args...; kwargs...)  kalman(Discrete, args...; kwargs...)

"""
    place(A, B, p, opt=:c; direct = false)
    place(sys::StateSpace, p, opt=:c; direct = false)

Calculate the gain matrix `K` such that `A - BK` has eigenvalues `p`.

    place(A, C, p, opt=:o)
    place(sys::StateSpace, p, opt=:o)

Calculate the observer gain matrix `L` such that `A - LC` has eigenvalues `p`.

If `direct = true` and `opt = :o`, the the observer gain `K` is calculated such that `A - KCA` has eigenvalues `p`, this option is to be used together with `direct = true` in [`observer_controller`](@ref). 

Note: only apply `direct = true` to discrete-time systems.

Ref: "Computer-Controlled Systems" pp 140.

Uses Ackermann's formula for SISO systems and [`place_knvd`](@ref) for MIMO systems. 

Please note that this function can be numerically sensitive, solving the placement problem in extended precision might be beneficial.
"""
function place(A, B, p, opt=:c; direct = false, kwargs...)
    n = length(p)
    n != size(A,1) && error("Must specify as many poles as states")
    if opt === :c
        direct && error("direct = true only applies to observer design")
        n != size(B,1) && error("A and B must have same number of rows")
        if size(B,2) == 1
            acker(A, B, p)
        else
            place_knvd(A, B, p; kwargs...)
        end
    elseif opt === :o
        C = B # B is really the "C matrix"
        if direct
            C = C*A
        end
        n != size(C,2) && error("A and C must have same number of columns")
        if size(C,1) == 1
            acker(A', C', p)'
        else
            place_knvd(A', C', p; kwargs...)'
        end
    else
        error("fourth argument must be :c or :o")
    end
end
function place(sys::AbstractStateSpace, p, opt=:c; direct = false, kwargs...)
    if opt === :c
        return place(sys.A, sys.B, p, opt; kwargs...)
    elseif opt === :o
        iscontinuous(sys) && direct && error("direct = true only applies to discrete-time systems")
        return place(sys.A, sys.C, p, opt; direct, kwargs...)
    else
        error("third argument must be :c or :o")
    end
end


"""
    acker(A,B,P)

Implements Ackermann's formula for placing poles of (A-BK) in p

Ackermann's formula works for SISO systems only, but a trick is possible to make Ackermann work for MIMO systems:
The code below introduces a random projection matrix `P` that projects the input space to one dimension, and then shifts the application of `P` from `B` to `K`. 

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

See also [`place_knvd`](@ref) which naturally handles MIMO systems.
"""
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


"""
    place_knvd(A::AbstractMatrix, B, λ; verbose = false, init = :s)

Robust pole placement using the algorithm from
> "Robust Pole Assignment in Linear State Feedback", Kautsky, Nichols, Van Dooren

This implementation uses "method 0" for the X-step and the QR factorization for all factorizations.

This function will be called automatically when [`place`](@ref) is called with a MIMO system.

# Arguments:
- `init`: Determines the initialization strategy for the iterations for find the `X` matrix. Possible choices are `:id` (default), `:rand`, `:s`. 
"""
function place_knvd(A::AbstractMatrix, B, λ; verbose=false, init=:s, method = 0)
    n, m = size(B)
    T = float(promote_type(eltype(A), eltype(B)))
    CT = Complex{real(T)}
    λ = sort(λ, by=LinearAlgebra.eigsortby)
    length(λ) == size(A, 1) == n || error("Must specify as many poles as states")
    Λ = diagm(λ)
    QRB = qr(B)
    U0, U1 = QRB.Q[:, 1:m], QRB.Q[:, m+1:end] # TODO: check dimension
    Z = QRB.R
    R = svdvals(B)
    m = count(>(100*eps()*R[1]), R) # Rank of B
    if m == n # Easy case, B is full rank
        r = count(e->imag(e) == 0, λ)
        ABF = diagm(real(λ))
        j = r+1
        while j <= n-1
            ABF[j, j+1] = imag(λ[j])
            ABF[j+1, j] = - imag(λ[j])
            j += 2
        end;
        return B\(A - ABF) # Solve for F in (A - BF) = Λ
    end

    S = Matrix{CT}[]
    for j = 1:n
        qj = qr((U1'*(A- λ[j]*I))')
        # Ŝj = qj.Q[:, 1:n-m] # Needed for method 2
        Sj = qj.Q[:, n-m+1:n]
        push!(S, Sj)
    end
    if m == 1 # Shortcut
        verbose && @info "Shortcut"
        X = reduce(hcat, S)
        M = X*Λ/X
        F = real((Z\U0')*(M-A))
        return -F
    end

    # Init
    if init === :id
        X = Matrix(one(CT)*I(n))
    elseif init === :rand
        X = randn(CT, n, n)
    elseif init === :s
        X = zeros(CT, n, n)
        @views for j = 1:n
            X[:,j] = sum(S[j], dims=2)
            X[:,j] ./= norm(X[:,j])
        end
    else
        error("Unknown init method")
    end

    cond_old = float(T)(Inf)
    if method == 0
        for i = 1:200
            verbose && @info "Iteration $i"
            for j = 1:n
                Xj = qr(X[:, setdiff(1:n, j)])
                ỹ = Xj.Q[:, end]
                STy = S[j]'ỹ
                xj = S[j]*(STy ./ norm(STy))
                any(!isfinite, xj) && error("Not finite")
                X[:, j] = xj
            end
            c = cond(X)
            verbose && @info "cond(X) = $c"
            if cond_old - c < 1e-14
                break
            end
            cond_old = c
            i == 200 && @warn "Max iterations reached"
        end
    # elseif method == 1
    #     for i = 1:200
    #         verbose && @info "Iteration $i"
    #         for j = 1:n
    #             Xj = qr(X[:, setdiff(1:n, j)])
    #             Rj = Xj.R
    #             Qj = Xj.Q[:, 1:end-1]
    #             qj = Xj.Q[:, end]
    #             σpt = qj'S[j]
    #             σ = norm(σpt)
    #             p = vec(σpt) ./ σ
    #             ρ = sqrt(σ^(-2)*(w̃'w̃ + 1))
    #             xj = inv(ρ*σ)*S[j]

    #             any(!isfinite, xj) && error("Not finite")
    #             X[:, j] = xj
    #         end
    #         c = cond(X)
    #         verbose && @info "cond(X) = $c"
    #         if cond_old - c < 1e-14
    #             break
    #         end
    #         cond_old = c
    #         i == 200 && @warn "Max iterations reached"
    #     end
    else
        error("Only method 0 is implemented")
    end
    verbose && @info "norm X = $(norm(X))"
    M = X*Λ/X
    F = real((Z\U0')*(M-A))
    -F # Paper assumes positive feedback
end
