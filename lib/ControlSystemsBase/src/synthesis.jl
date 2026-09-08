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
    alpha_beta(alpha, Ts; beta = 2(2 - alpha) - 4√(1 - alpha))

Make an α-β tracker. The tracker estimates the position and the rate of a target. The function
returns a discrete-time `StateSpace` system. The input of the system is the measured position.
The outputs are the estimates ``[x̂, v̂]``.

```math
\\begin{aligned}
r(k) &= y(k) - \\big(x̂(k-1) + T_s v̂(k-1)\\big) \\\\
x̂(k) &= x̂(k-1) + T_s v̂(k-1) + α\\, r(k) \\\\
v̂(k) &= v̂(k-1) + \\dfrac{β}{T_s} r(k)
\\end{aligned}
```

The system is [`observer_filter`](@ref) of a double integrator with the gain
``K = [α,\\ β/T_s]``. The state of the system is the a-posteriori estimate. Thus the system is
strictly proper. The rate output is a filtered derivative of the input. A separate differentiator
is not necessary.

!!! note "The sample index"
    The state is ``x̂(k|k)``. This estimate includes the measurement ``y(k)``. The state-space
    update is ``x(k+1) = Ax(k) + Bu(k)``. Therefore, in a simulation, this estimate is at output
    index ``k+1``. The output at index ``k`` includes the measurements up to ``y(k-1)``.
    [`observer_filter`](@ref) uses the same convention.

`alpha` sets the effect of each measurement on the position estimate. The value of `alpha` must
be in the range ``0 < α < 1``. The value of `beta` must be in the range ``0 < β ≤ 2 - α``. The
default `beta` is Kalata's steady-state relation ``β = 2(2 - α) - 4\\sqrt{1 - α}``. This relation
makes the filter equal to the steady-state Kalman filter for a target that has a constant rate.
Thus `alpha` is a sufficient tuning parameter.

For a target that accelerates, use [`alpha_beta_gamma`](@ref). If you know the noise covariances
and not the gains, use [`kalman`](@ref) with [`observer_filter`](@ref).

# Example
The rate estimate converges to the slope of a ramp input.
```jldoctest
julia> using ControlSystemsBase

julia> Ts = 0.1; sys = alpha_beta(0.5, Ts);

julia> size(sys)
(2, 1)

julia> res = lsim(sys, (x, t) -> [t], 0:Ts:5);

julia> round(res.y[2, end], digits = 3)
1.0
```

# Extended help
The estimation-error dynamics are ``(I - KC)A``. The gains set the eigenvalues of this matrix.
The default `beta` gives a complex pole pair for each value of `alpha`. To put both error poles
on the real axis at the same radius ``s ∈ (0, 1)``, give both gains:

``α = 1 - s^2``, ``β = (1 - s)^2``.

These gains make the filter critically damped. The error then decays with no oscillation.

Do not give only `alpha`. If you give only `alpha`, `beta` keeps its default value, and the poles
do not stay on the real axis.

``s`` replaces `alpha` as the tuning parameter. The error decays as ``k s^k``. Use
``-T_s/\\ln s`` as an estimate of the time constant. A smaller value of ``s`` gives faster
tracking and more measurement noise in the estimates.
"""
function alpha_beta(alpha, Ts; beta = 2 * (2 - alpha) - 4 * sqrt(1 - alpha))
    0 < alpha < 1 || throw(ArgumentError("alpha must satisfy 0 < alpha < 1, got $alpha"))
    Ts > 0 || throw(ArgumentError("Ts must be positive, got $Ts"))
    T = float(promote_type(typeof(alpha), typeof(beta), typeof(Ts)))
    A = T[1 Ts; 0 1]
    C = T[1 0]
    K = T[alpha; beta/Ts;;]
    ss((I - K * C) * A, K, Matrix{T}(I, 2, 2), zeros(T, 2, 1), Ts)
end

"""
    alpha_beta_gamma(alpha, Ts; beta = 2(2 - alpha) - 4√(1 - alpha), gamma = beta^2 / (2alpha))

Make an α-β-γ tracker. The tracker estimates the position, the rate and the acceleration of a
target. The function returns a discrete-time `StateSpace` system. The input of the system is the
measured position. The outputs are the estimates ``[x̂, v̂, â]``.

```math
\\begin{aligned}
r(k) &= y(k) - \\big(x̂(k-1) + T_s v̂(k-1) + \\tfrac{T_s^2}{2} â(k-1)\\big) \\\\
x̂(k) &= x̂(k-1) + T_s v̂(k-1) + \\tfrac{T_s^2}{2} â(k-1) + α\\, r(k) \\\\
v̂(k) &= v̂(k-1) + T_s â(k-1) + \\dfrac{β}{T_s} r(k) \\\\
â(k) &= â(k-1) + \\dfrac{γ}{T_s^2} r(k)
\\end{aligned}
```

The system is [`observer_filter`](@ref) of a triple integrator with the gain
``K = [α,\\ β/T_s,\\ γ/T_s^2]``. The rule for the sample index that applies to
[`alpha_beta`](@ref) also applies to this filter.

An α-β tracker predicts with a constant rate. If the target accelerates, the rate estimate of an
α-β tracker has an error. This error is proportional to the acceleration. No value of `alpha` and
`beta` removes this error. This filter also predicts with the acceleration. Thus the error is not
present.

Use this filter when the sample rate is much higher than the frequency content of the signal. In
this condition, a long filter memory causes only a small lag. The filter then decreases the
effect of the measurement noise.

The default gains are Kalata's steady-state relations. These gains make the filter equal to the
steady-state Kalman filter for a target that has a constant acceleration. Thus `alpha` is a
sufficient tuning parameter. If ``α = 0.5``, then ``β ≈ 0.1716`` and ``γ ≈ 0.0294``.

# Example
```jldoctest
julia> using ControlSystemsBase

julia> sys = alpha_beta_gamma(0.5, 0.1);

julia> size(sys)
(3, 1)
```

# Extended help
The default gains are not critically damped. This is also true for [`alpha_beta`](@ref). The
Kalata gains give a complex pole pair for each value of `alpha`. To put all three error poles on
the real axis at the same radius ``s ∈ (0, 1)``, give all three gains:

``α = 1 - s^3``, ``β = \\tfrac{3}{2}(1 - s)^2(1 + s)``, ``γ = (1 - s)^3``.

Do not give only `alpha`. If you give only `alpha`, `beta` and `gamma` keep their default values,
and the poles do not stay on the real axis. For example, at ``s = 0.9`` the critically damped
gains are ``(0.271, 0.0285, 0.001)``. If you give `alpha = 0.271` and no other gain, you get
``β = 0.0427`` and ``γ = 0.0034``.

The pole has a multiplicity of three. Thus the error decays as ``k^2 s^k``, and the settling time
is longer than the pole radius alone indicates. Use ``-T_s/\\ln s`` as an estimate of the time
constant. At ``s = 0.9`` this estimate is approximately 9.5 samples.

These relations apply to the gain convention that is shown above. In this convention, the
acceleration correction is ``γ/T_s^2`` and the rate correction is ``β/T_s``. Other references put
``2γ`` in this position. The relations in those references are thus different by a factor of two.
"""
function alpha_beta_gamma(alpha, Ts; beta = 2 * (2 - alpha) - 4 * sqrt(1 - alpha),
                                     gamma = beta^2 / (2 * alpha))
    0 < alpha < 1 || throw(ArgumentError("alpha must satisfy 0 < alpha < 1, got $alpha"))
    Ts > 0 || throw(ArgumentError("Ts must be positive, got $Ts"))
    T = float(promote_type(typeof(alpha), typeof(beta), typeof(gamma), typeof(Ts)))
    A = T[1 Ts Ts^2/2; 0 1 Ts; 0 0 1]
    C = T[1 0 0]
    K = T[alpha; beta/Ts; gamma/Ts^2;;]
    ss((I - K * C) * A, K, Matrix{T}(I, 3, 3), zeros(T, 3, 1), Ts)
end

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
    n != size(A,1) && error("Must specify as many poles as the state dimension")
    L = if opt === :c
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
    if isreal(A) && isreal(B) && is_self_conjugate(p)
        @assert all(abs(imag(l)) .< 1e-6 for l in L) "Expected real coefficient in feedback gain, got complex: $L"
        return real(L)
    end
    return L
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
- `init`: Determines the initialization strategy for the iterations for find the `X` matrix. Possible choices are `:id`, `:rand`, `:s` (default). 
"""
function place_knvd(A::AbstractMatrix, B, λ; verbose=false, init=:s, method = 0)
    n, m = size(B)
    T = float(promote_type(eltype(A), eltype(B)))
    CT = Complex{real(T)}
    λ = sort(vec(λ), by=LinearAlgebra.eigsortby)
    length(λ) == size(A, 1) == n || error("Must specify as many poles as the state dimension")
    Λ = diagm(λ)
    R = svdvals(B)
    m = count(>(100*eps()*R[1]), R) # Rank of B
    QRB = qr(B, ColumnNorm())
    U0, U1 = QRB.Q[:, 1:m], QRB.Q[:, m+1:end] # TODO: check dimension
    Z = (QRB.R*QRB.P')[:, 1:m] 
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

    mB = size(B, 2)
    if mB > m
        # several inputs but not full column rank, this case must be handled separately
        # when B does not have full column rank but that rank is not 1. In that case, find B2 and T from rank-revealing QR (qr(B, ColumnNorm())
        verbose && @info "Projecting down to rank of B"
        B2 = QRB.Q[:, 1:m]
        T = QRB.P * QRB.R[1:m, :]'
        F = place(A, B2, λ; verbose, init, method)
        return pinv(T)'*F
    end

    S = Matrix{CT}[]
    for j = 1:n
        H = (U1'*(A- λ[j]*I))
        qj = qr(H')
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
