"""
    are(::Continuous, A, B, Q, R)

Compute 'X', the solution to the continuous-time algebraic Riccati equation,
defined as A'X + XA - (XB)R^-1(B'X) + Q = 0, where R is non-singular.

Uses `MatrixEquations.arec`.
"""
function are(::ContinuousType, A::AbstractMatrix, B, Q, R)
    arec(A, B, R, Q)[1]
end

"""
    are(::Discrete, A, B, Q, R; kwargs...)

Compute `X`, the solution to the discrete-time algebraic Riccati equation,
defined as A'XA - X - (A'XB)(B'XB + R)^-1(B'XA) + Q = 0, where Q>=0 and R>0

Uses `MatrixEquations.ared`. For keyword arguments, see the docstring of `ControlSystems.MatrixEquations.ared`
"""
function are(::DiscreteType, A::AbstractMatrix, B, Q, R; kwargs...)
    ared(A, B, R, Q; kwargs...)[1]
end

are(t::TimeEvolType, A::Number, B::Number, Q::Number, R::Number) = are(t, fill(A,1,1),fill(B,1,1),fill(Q,1,1),fill(R,1,1))

@deprecate care(args...; kwargs...) are(Continuous, args...; kwargs...)
@deprecate dare(args...; kwargs...) are(Discrete, args...; kwargs...)

"""
    lyap(A, Q; kwargs...)

Compute the solution `X` to the discrete Lyapunov equation
`AXA' - X + Q = 0`.

Uses `MatrixEquations.lyapc / MatrixEquations.lyapd`. For keyword arguments, see the docstring of `ControlSystems.MatrixEquations.lyapc / ControlSystems.MatrixEquations.lyapd`
"""
function LinearAlgebra.lyap(::DiscreteType, A::AbstractMatrix, Q; kwargs...)
    lyapd(A, Q; kwargs...)
end

LinearAlgebra.lyap(::ContinuousType, args...; kwargs...) = lyapc(args...; kwargs...)
LinearAlgebra.lyap(::DiscreteType, args...; kwargs...) = lyapd(args...; kwargs...)

plyap(::ContinuousType, args...; kwargs...) = MatrixEquations.plyapc(args...; kwargs...)
plyap(::DiscreteType, args...; kwargs...) = MatrixEquations.plyapd(args...; kwargs...)

@deprecate dlyap(args...; kwargs...) lyap(Discrete, args...; kwargs...)


"""
    U = grampd(sys, opt; kwargs...)

Return a Cholesky factor `U` of the grammian of system `sys`. If `opt` is `:c`, computes the
controllability grammian `G = U*U'`. If `opt` is `:o`, computes the observability
grammian `G = U'U`.

Obtain a `Cholesky` object by `Cholesky(U)` for observability grammian

Uses `MatrixEquations.plyapc/plyapd`. For keyword arguments, see the docstring of `ControlSystems.MatrixEquations.plyapc/plyapd`
"""
function grampd(sys::AbstractStateSpace, opt::Symbol; kwargs...)
    if !isstable(sys)
        error("gram only valid for stable A")
    end
    if opt === :c
        plyap(sys.timeevol, sys.A, sys.B; kwargs...)
    elseif opt === :o
        plyap(sys.timeevol, sys.A', sys.C'; kwargs...)
    else
        error("opt must be either :c for controllability grammian, or :o for
                observability grammian")
    end
end

"""
    gram(sys, opt; kwargs...)

Compute the grammian of system `sys`. If `opt` is `:c`, computes the
controllability grammian. If `opt` is `:o`, computes the observability
grammian.

See also [`grampd`](@ref)
For keyword arguments, see [`grampd`](@ref).
"""
function gram(sys::AbstractStateSpace, opt::Symbol; kwargs...)
    U = grampd(sys, opt; kwargs...)
    opt === :c ? U*U' : U'U
end

"""
    obsv(A, C, n=size(A,1))
    obsv(sys, n=sys.nx)

Compute the observability matrix with `n` rows for the system described by `(A, C)` or `sys`. Providing the optional `n > sys.nx` returns an extended observability matrix.

Note that checking for observability by computing the rank from `obsv` is
not the most numerically accurate way, a better method is checking if
`gram(sys, :o)` is positive definite.
"""
function obsv(A::AbstractMatrix, C::AbstractMatrix, n::Int = size(A,1))
    T = promote_type(eltype(A), eltype(C))
    nx = size(A, 1)
    ny = size(C, 1)
    if nx != size(C, 2)
        throw(ArgumentError("C must have the same number of columns as A"))
    end
    res = fill(zero(T), n*ny, nx)
    res[1:ny, :] = C
    for i=1:n-1
        res[(1 + i*ny):(1 + i)*ny, :] = res[((i - 1)*ny + 1):i*ny, :] * A
    end
    return res
end
obsv(sys::AbstractStateSpace, n::Int = sys.nx) = obsv(sys.A, sys.C, n)

"""`ctrb(A, B)` or `ctrb(sys)`

Compute the controllability matrix for the system described by `(A, B)` or
`sys`.

Note that checking for controllability by computing the rank from
`ctrb` is not the most numerically accurate way, a better method is
checking if `gram(sys, :c)` is positive definite."""
function ctrb(A::AbstractMatrix, B::AbstractMatrix)
    T = promote_type(eltype(A), eltype(B))
    n = size(A, 1)
    nu = size(B, 2)
    if n != size(B, 1)
        throw(ArgumentError("B must have the same number of rows as A"))
    end
    res = fill(zero(T), n, n*nu)
    res[:, 1:nu] = B
    for i=1:n-1
        res[:, (1 + i*nu):(1 + i)*nu] = A * res[:, ((i - 1)*nu + 1):i*nu]
    end
    return res
end
ctrb(sys::AbstractStateSpace) = ctrb(sys.A, sys.B)

"""`P = covar(sys, W)`

Calculate the stationary covariance `P = E[y(t)y(t)']` of the output `y` of a
`StateSpace` model `sys` driven by white Gaussian noise `w` with covariance
`E[w(t)w(τ)]=W*δ(t-τ)` (δ is the Dirac delta).

Remark: If `sys` is unstable then the resulting covariance is a matrix of `Inf`s.
Entries corresponding to direct feedthrough (D*W*D' .!= 0) will equal `Inf`
for continuous-time systems."""
function covar(sys::AbstractStateSpace, W)
    (A, B, C, D) = ssdata(sys)
    if !isa(W, UniformScaling) && (size(B,2) != size(W, 1) || size(W, 1) != size(W, 2))
        error("W must be a square matrix the same size as `sys.B` columns")
    end
    isa(W, UniformScaling) && (W = I(size(B, 2)))
    if !isstable(sys)
        return fill(Inf,(size(C,1),size(C,1)))
    end
    Wc = cholesky(W).L
    Q1 = sys.nx == 0 ? B*Wc : try
        plyap(sys.timeevol, A, B*Wc)
    catch e
        @error("No solution to the Lyapunov equation was found in covar.")
        rethrow(e)
    end
    P1 = C*Q1
    P = P1*P1'
    if iscontinuous(sys)
        #Variance and covariance infinite for direct terms
        direct_noise = D*W*D'
        for i in 1:size(C,1)
            if direct_noise[i,i] != 0
                P[i,:] .= Inf
                P[:,i] .= Inf
            end
        end
    else
        P += D*W*D'
    end
    return P
end

covar(sys::TransferFunction, W) = covar(ss(sys), W)


# Note: the H∞ norm computation is probably not as accurate as with SLICOT,
# but this seems to be still reasonably ok as a first step
"""
    norm(sys, p=2; tol=1e-6)

`norm(sys)` or `norm(sys,2)` computes the H2 norm of the LTI system `sys`.

`norm(sys, Inf)` computes the H∞ norm of the LTI system `sys`.
The H∞ norm is the same as the H∞ for stable systems, and Inf for unstable systems.
If the peak gain frequency is required as well, use the function `hinfnorm` instead.
See [`hinfnorm`](@ref) for further documentation.

`tol` is an optional keyword argument, used only for the computation of L∞ norms.
It represents the desired relative accuracy for the computed L∞ norm
(this is not an absolute certificate however).

`sys` is first converted to a `StateSpace` model if needed.
"""
function LinearAlgebra.norm(sys::AbstractStateSpace, p::Real=2; tol=1e-6)
    if p == 2
        return sqrt(max(0,tr(covar(sys, I))))
    elseif p == Inf
        return hinfnorm(sys; tol=tol)[1]
    else
        error("`p` must be either `2` or `Inf`")
    end
end
LinearAlgebra.norm(sys::TransferFunction, p::Real=2; tol=1e-6) = norm(ss(sys), p, tol=tol)


"""
    Ninf, ω_peak = hinfnorm(sys; tol=1e-6)

Compute the H∞ norm `Ninf` of the LTI system `sys`, together with a frequency
`ω_peak` at which the gain Ninf is achieved.

`Ninf := sup_ω σ_max[sys(iω)]`  if `G` is stable (σ_max = largest singular value)
      :=        `Inf'           if `G` is unstable

`tol` is an optional keyword argument for the desired relative accuracy for
the computed H∞ norm (not an absolute certificate).

`sys` is first converted to a state space model if needed.

The continuous-time L∞ norm computation implements the 'two-step algorithm' in:\\
**N.A. Bruinsma and M. Steinbuch**, 'A fast algorithm to compute the H∞-norm of
a transfer function matrix', Systems and Control Letters (1990), pp. 287-293.

For the discrete-time version, see:\\
**P. Bongers, O. Bosgra, M. Steinbuch**, 'L∞-norm calculation for generalized
state space systems in continuous and discrete time', American Control Conference, 1991.

See also [`linfnorm`](@ref).
"""
hinfnorm(sys::AbstractStateSpace{<:Continuous}; tol=1e-6) = _infnorm_two_steps_ct(sys, :hinf, tol)
hinfnorm(sys::AbstractStateSpace{<:Discrete}; tol=1e-6) = _infnorm_two_steps_dt(sys, :hinf, tol)
hinfnorm(sys::TransferFunction; tol=1e-6) = hinfnorm(ss(sys); tol=tol)

"""
    Ninf, ω_peak = linfnorm(sys; tol=1e-6)

Compute the L∞ norm `Ninf` of the LTI system `sys`, together with a frequency
`ω_peak` at which the gain `Ninf` is achieved.

`Ninf := sup_ω σ_max[sys(iω)]` (σ_max denotes the largest singular value)

`tol` is an optional keyword argument representing the desired relative accuracy for
the computed L∞ norm (this is not an absolute certificate however).

`sys` is first converted to a state space model if needed.

The continuous-time L∞ norm computation implements the 'two-step algorithm' in:\\
**N.A. Bruinsma and M. Steinbuch**, 'A fast algorithm to compute the H∞-norm of
a transfer function matrix', Systems and Control Letters (1990), pp. 287-293.

For the discrete-time version, see:\\
**P. Bongers, O. Bosgra, M. Steinbuch**, 'L∞-norm calculation for generalized
state space systems in continuous and discrete time', American Control Conference, 1991.

See also [`hinfnorm`](@ref).
"""
function linfnorm(sys::AbstractStateSpace; tol=1e-6)
    if iscontinuous(sys)
        return _infnorm_two_steps_ct(sys, :linf, tol)
    else
        return _infnorm_two_steps_dt(sys, :linf, tol)
    end
end
linfnorm(sys::TransferFunction; tol=1e-6) = linfnorm(ss(sys); tol=tol)

function _infnorm_two_steps_ct(sys::AbstractStateSpace, normtype::Symbol, tol=1e-6, maxIters=250, approximag=1e-10)
    # norm type :hinf or :linf the reason that to not use `hinfnorm(sys) = isstable(sys) : linfnorm ? (Inf, Nan)`
    # is to avoid re computing the poles and return the peak frequencies for, e.g., 1/(s^2 + 1)
    # `maxIters`: the maximum  number of iterations allowed in the algorithm (default 1000)
    # approximag is a tuning parameter: what does it mean for a number to be on the imaginary axis
    # Because of this tuning for example, the relative precision that we provide on the norm computation
    # is not a true guarantee, more an order of magnitude
    # outputs: An approximatation of the L∞ norm and the frequency ω_peak at which it is achieved
    # QUESTION: The tolerance for determining if there are poles on the imaginary axis
    # would not be very appropriate for systems with slow dynamics?
    T = promote_type(real(numeric_type(sys)), Float64)
    on_imag_axis = z -> abs(real(z)) < approximag # Helper fcn for readability

    if sys.nx == 0  # static gain
        return (T(opnorm(sys.D)), T(0))
    end

    pole_vec = poles(sys)

    # Check if there is a pole on the imaginary axis
    pidx = findfirst(on_imag_axis, pole_vec)
    if !(pidx isa Nothing)
        return (T(Inf), T(imag(pole_vec[pidx])))
        # note: in case of cancellation, for s/s for example, we return Inf, whereas Matlab returns 1
    end

    if normtype === :hinf && any(z -> real(z) > 0, pole_vec)
        return T(Inf), T(NaN) # The system is unstable
    end

    # Initialization: computation of a lower bound from 3 terms
    if isreal(pole_vec)  # only real poles
        ω_p = minimum(abs.(pole_vec))
    else  # at least one pair of complex poles
        maxidx = argmax([abs(imag(p)/real(p))/abs(p) for p in pole_vec])
        ω_p = abs(pole_vec[maxidx])
    end

    m_vec_init = [0, ω_p, Inf]

    (lb, idx) = findmax([opnorm(evalfr(sys, im*m_vec_init[1]));
                         opnorm(evalfr(sys, im*m_vec_init[2]));
                         opnorm(sys.D)])
    ω_peak = m_vec_init[idx]
    lb == 0 && (return zero(T), zero(T))
    # Iterations
    for iter=1:maxIters
        gamma = (1+2*T(tol))*lb
        R = sys.D'*sys.D - gamma^2*I
        S = sys.D*sys.D' - gamma^2*I
        M = sys.A-sys.B*(R\sys.D')*sys.C
        H = [         M              -gamma*sys.B*(R\sys.B') ;
               gamma*sys.C'*(S\sys.C)            -M'            ]

        Λ = complex(eigvals(H)) # To make type stable

        if numeric_type(sys) <: Real
            # Only need to consider one eigenvalue in each complex-conjugate pairs
            filter!(z -> imag(z) >= 0, Λ)
        end

        # Find eigenvalues on the imaginary axis
        Λ_on_imag_axis = filter(on_imag_axis, Λ)

        ω_vec = imag.(Λ_on_imag_axis)

        sort!(ω_vec)

        if isempty(ω_vec)
            return T((1+tol)*lb), T(ω_peak)
        end

        # Improve the lower bound
        # if not empty, ω_vec contains at least two values
        for k=1:length(ω_vec)-1
            mk = (ω_vec[k] + ω_vec[k+1])/2
            sigmamax_mk = opnorm(evalfr(sys,mk*1im))
            if sigmamax_mk > lb
                lb = sigmamax_mk
                ω_peak = mk
            end
        end
    end
    error("In _infnorm_two_steps_dt: The computation of the H∞/L∞ norm did not converge in $maxIters iterations")
end

function _infnorm_two_steps_dt(sys::AbstractStateSpace, normtype::Symbol, tol=1e-6, maxIters=250, approxcirc=1e-8)
    # Discrete-time version of linfnorm_two_steps_ct above
    # Compuations are done in normalized frequency θ

    on_unit_circle = z -> abs(abs(z) - 1) < approxcirc # Helper fcn for readability

    T = promote_type(real(numeric_type(sys)), Float64, typeof(true/sys.Ts))
    Tw = typeof(one(T)/sys.Ts)

    if sys.nx == 0  # static gain
        return (T(opnorm(sys.D)), Tw(0))
    end

    pole_vec = poles(sys)

    # Check if there is a pole on the unit circle
    pidx = findfirst(on_unit_circle, pole_vec)
    if !(pidx isa Nothing)
        return T(Inf), Tw(angle(pole_vec[pidx])/sys.Ts)
    end

    if normtype == :hinf && any(z -> abs(z) > 1, pole_vec)
        return T(Inf), Tw(NaN) # The system is unstable
    end

    # Initialization: computation of a lower bound from 3 terms

    if isreal(pole_vec)  # not just real poles
        # find frequency of pôle closest to unit circle
        θ_p = angle(pole_vec[argmin(abs.(abs.(pole_vec).-1))])
    else
        θ_p = T(pi)/2
    end

    if isreal(pole_vec)  # only real poles
        ω_p = minimum(abs.(pole_vec))
    else  # at least one pair of complex poles
        maxidx = argmax([abs(imag(p)/real(p))/abs(p) for p in pole_vec])
        ω_p = abs(pole_vec[maxidx])
    end

    m_vec_init = [0, θ_p, pi]

    (lb, idx) = findmax([opnorm(evalfr(sys, exp(im*m))) for m in m_vec_init])
    θ_peak = m_vec_init[idx]

    # Iterations
    for iter=1:maxIters
        gamma = (1+2*T(tol))*lb
        R = gamma^2*I - sys.D'*sys.D
        RinvDt = R\sys.D'
        L = [ sys.A+sys.B*RinvDt*sys.C  sys.B*(R\sys.B');
              zeros(T, sys.nx,sys.nx)      I]
        M = [ I                                 zeros(T, sys.nx,sys.nx);
              sys.C'*(I+sys.D*RinvDt)*sys.C     L[1:sys.nx,1:sys.nx]']

        Λ = complex(eigvals(L,M)) # complex is to ensure type stability

        if numeric_type(sys) <: Real
            # Only need to consider one eigenvalue in each complex-conjugate pairs
            filter!(z -> imag(z) >= 0, Λ)
        end

        # Find eigenvalues on the unit circle
        Λ_on_unit_cirlce = filter(on_unit_circle, Λ)

        θ_vec = angle.(Λ_on_unit_cirlce)

        sort!(θ_vec)

        if isempty(θ_vec)
            return T((1+tol)*lb), Tw(θ_peak/sys.Ts)
        end

        # Improve the lower bound
        # if not empty, θ_vec contains at least two values
        for k=1:length(θ_vec)-1
            mk = (θ_vec[k] + θ_vec[k+1])/2
            sigmamax_mk = opnorm(evalfr(sys,exp(mk*1im)))
            if sigmamax_mk > lb
                lb = sigmamax_mk
                θ_peak = mk
            end
        end
    end
    error("In _infnorm_two_steps_dt: The computation of the L∞ norm did not converge in $maxIters iterations")
end


"""
    S, P, B = balance(A[, perm=true])

Compute a similarity transform `T = S*P` resulting in `B = T\\A*T` such that the row
and column norms of `B` are approximately equivalent. If `perm=false`, the
transformation will only scale `A` using diagonal `S`, and not permute `A` (i.e., set `P=I`).
"""
function balance(A, perm::Bool=true)
    n = LinearAlgebra.checksquare(A)
    B = copy(A)
    job = perm ? 'B' : 'S'
    ilo, ihi, scaling = LAPACK.gebal!(job, B)

    S = diagm(0 => scaling)
    for j = 1:(ilo-1)   S[j,j] = 1 end
    for j = (ihi+1):n   S[j,j] = 1 end

    P = Matrix{Int}(I,n,n)
    if perm
        if ilo > 1
            for j = (ilo-1):-1:1 cswap!(j, round(Int, scaling[j]), P) end
        end
        if ihi < n
            for j = (ihi+1):n    cswap!(j, round(Int, scaling[j]), P) end
        end
    end
    return S, P, B
end

function cswap!(i::Integer, j::Integer, X::StridedMatrix)
    for k = 1:size(X,1)
        X[i, k], X[j, k] = X[j, k], X[i, k]
    end
end



"""
`sysr, G, T = balreal(sys::StateSpace)`

Calculates a balanced realization of the system sys, such that the observability and reachability gramians of the balanced system are equal and diagonal `G`. `T` is the similarity transform between the old state `x` and the new state `z` such that `Tz = x`.

See also `gram`, `baltrunc`

Reference: Varga A., Balancing-free square-root algorithm for computing singular perturbation approximations.
"""
function balreal(sys::ST) where ST <: AbstractStateSpace
    # This code is adapted from DescriptorSystems.jl
    # originally written by Andreas Varga
    # https://github.com/andreasvarga/DescriptorSystems.jl/blob/dd144828c3615bea2d5b4977d7fc7f9677dfc9f8/src/order_reduction.jl#L622
    # with license https://github.com/andreasvarga/DescriptorSystems.jl/blob/main/LICENSE.md
    A,B,C,D = ssdata(sys)

    SF = schur(A)
    bs = SF.Z'*B
    cs = C*SF.Z

    S = MatrixEquations.plyaps(SF.T, bs; disc = isdiscrete(sys))
    R = MatrixEquations.plyaps(SF.T', cs'; disc = isdiscrete(sys))
    SV = svd!(R*S)

    U,Σ,V = SV

    # Determine the order of a minimal realization to √ϵ tolerance
    rmin = count(Σ .> sqrt(eps())*Σ[1])
    i1 = 1:rmin
    Σ = Σ[i1]


    hsi2 = Diagonal(1 ./sqrt.(Σ))
    L = lmul!(R',view(U,:,i1))*hsi2
    Tr = lmul!(S,V[:,i1])*hsi2
    # return the minimal balanced system
    T = L'SF.Z'
    return ss(L'SF.T*Tr, L'bs, cs*Tr, sys.D, sys.timeevol), Diagonal(Σ), T
end


"""
    sysr, G, T = baltrunc(sys::StateSpace; atol = √ϵ, rtol=1e-3, n = nothing, residual = false)

Reduces the state dimension by calculating a balanced realization of the system sys, such that the observability and reachability gramians of the balanced system are equal and diagonal `G`, and truncating it to order `n`. If `n` is not provided, it's chosen such that all states corresponding to singular values less than `atol` and less that `rtol σmax` are removed.

`T` is the similarity transform between the old state `x` and the newstate `z` such that `Tz = x`.

If `residual = true`, matched static gain is achieved through "residualization", i.e., setting
```math
0 = A_{21}x_{1} + A_{22}x_{2} + B_{2}u
```
where indices 1/2 correspond to the remaining/truncated states respectively.

See also `gram`, `balreal`

Glad, Ljung, Reglerteori: Flervariabla och Olinjära metoder
"""
function baltrunc(sys::ST; atol = sqrt(eps()), rtol = 1e-3, n = nothing, residual=false) where ST <: AbstractStateSpace
    sysbal, S, T = balreal(sys)
    S = diag(S)
    if n === nothing
        S = S[S .>= atol]
        S = S[S .>= S[1]*rtol]
        n = length(S)
    else
        S = S[1:n]
    end
    i1 = 1:n
    if residual
        A,B,C,D = ssdata(sysbal)
        i2 = n+1:size(A, 1)
        A11 = A[i1, i1]
        A12 = A[i1, i2]
        A21 = A[i2, i1]
        A22 = -A[i2, i2]
        isdiscrete(sys) && (A22 += I)
        B1 = B[i1, :]
        B2 = B[i2, :]
        C1 = C[:, i1]
        C2 = C[:, i2]
        A2221 = A22\A21
        A = A11 + A12*(A2221)
        B = B1 + (A12/A22)*B2
        C = C1 + C2*A2221
        D = D + (C2/A22)*B2
    else 
        A = sysbal.A[i1,i1]
        B = sysbal.B[i1,:]
        C = sysbal.C[:,i1]
        D = sysbal.D
    end

    return ST(A,B,C,D,sys.timeevol), diagm(S), T
end

"""
    syst = similarity_transform(sys, T; unitary=false)
Perform a similarity transform `T : Tx̃ = x` on `sys` such that
```
Ã = T⁻¹AT
B̃ = T⁻¹ B
C̃ = CT
D̃ = D
```

If `unitary=true`, `T` is assumed unitary and the matrix adjoint is used instead of the inverse.
See also [`balance_statespace`](@ref).
"""
function similarity_transform(sys::ST, T; unitary=false) where ST <: AbstractStateSpace
    if unitary
        A = T'sys.A*T
        B = T'sys.B
    else
        Tf = factorize(T)
        A = Tf\sys.A*T
        B = Tf\sys.B
    end
    C = sys.C*T
    D = sys.D
    ST(A,B,C,D,sys.timeevol)
end


"""
    sysi = innovation_form(sys, R1, R2[, R12])
    sysi = innovation_form(sys; sysw=I, syse=I, R1=I, R2=I)

Takes a system
```
x' = Ax + Bu + w ~ R1
y  = Cx + Du + e ~ R2
```
and returns the system
```
x' = Ax + Kv
y  = Cx + v
```
where `v` is the innovation sequence.

If `sysw` (`syse`) is given, the covariance resulting in filtering noise with `R1` (`R2`) through `sysw` (`syse`) is used as covariance.

See Stochastic Control, Chapter 4, Åström
"""
function innovation_form(sys::ST, R1, R2, args...) where ST <: AbstractStateSpace
    K = kalman(sys, R1, R2, args...)
    innovation_form(sys, K)
end
# Set D = I to get transfer function H = I + C(sI-A)\ K
function innovation_form(sys::ST; sysw=I, syse=I, R1=I, R2=I) where ST <: AbstractStateSpace
	K = kalman(sys, covar(sysw,R1), covar(syse, R2))
	ST(sys.A, K, sys.C, Matrix{eltype(sys.A)}(I, sys.ny, sys.ny), sys.timeevol)
end


"""
    sysi = innovation_form(sys, K)    

Takes a system
```
x' = Ax + Bu + Kv
y  = Cx + Du + v
```
and returns the system
```
x' = Ax + Kv
y  = Cx + v
```
where `v` is the innovation sequence.

See Stochastic Control, Chapter 4, Åström
"""
function innovation_form(sys::ST, K) where ST <: AbstractStateSpace
    ss(sys.A, K, sys.C, Matrix{eltype(sys.A)}(I, sys.ny, sys.ny), sys.timeevol)
end

"""
    observer_predictor(sys::AbstractStateSpace, K; h::Int = 1)
    observer_predictor(sys::AbstractStateSpace, R1, R2[, R21])

If `sys` is continuous, return the observer predictor system
x̂' = (A - KC)x̂ + (B-KD)u + Ky
ŷ  = Cx + Du
with the input equation [B-KD K] * [u; y]

If `sys` is discrete, the prediction horizon `h` may be specified, in which case measurements up to and including time `t-h` and inputs up to and including time `t` are used to predict `y(t)`.

If covariance matrices `R1, R2` are given, the kalman gain `K` is calculated using [`kalman`](@ref).

See also [`innovation_form`](@ref) and [`observer_controller`](@ref).
"""
function observer_predictor(sys::AbstractStateSpace, R1, R2::Union{AbstractArray, UniformScaling}, args...; kwargs...)
    K = kalman(sys, R1, R2, args...)
    observer_predictor(sys, K; kwargs...)
end

function observer_predictor(sys::AbstractStateSpace, K::AbstractMatrix; h::Integer = 1)
    h >= 1 || throw(ArgumentError("h must be positive."))
    ny = noutputs(sys)
    size(K, 1) == sys.nx && size(K,2) == ny || throw(ArgumentError("K has the wrong size, expected $((sys.nx, ny))"))
    A,B,C,D = ssdata(sys)
    if h == 1
        ss(A-K*C, [B-K*D K], C, [D zeros(ny, ny)], sys.timeevol)
    else
        isdiscrete(sys) || throw(ArgumentError("A prediction horizon is only supported for discrete systems. "))
        # The impulse response of the innovation form calculates the influence of a measurement at time t on the prediction at time t+h
        # Below, we form a system del (delay) that convolves the input (y) with the impulse response
        # We then add the output again to account for the fact that we propagated error and not measurement
        inn = innovation_form(sys, K)
        ts = (0:h-1) .* sys.Ts
        imp = impulse(inn, ts).y * sys.Ts # Impulse response differs from Markov params by 1/Ts
        del_components = map(Iterators.product(1:inn.ny, 1:inn.nu)) do (i,j)
            tf(imp[i,:,j], [1; zeros(h-1)]) # This forms a system that convolves the input with the impulse response
        end
        del = tf(first.(numvec.(del_components)), first.(denvec.(del_components)), sys.timeevol) |> ss
        pe = ss(A-K*C, [B-K*D K], C, [D -I(ny)], sys.timeevol) # prediction error system ŷ-y
        return ss([zero(D) I(ny)], sys.timeevol) + del*pe # add y back to compensate for -y in pe
    end
end

"""
    cont = observer_controller(sys, L::AbstractMatrix, K::AbstractMatrix)

Return the observer_controller `cont` that is given by
`ss(A - B*L - K*C + K*D*L, K, L, 0)`

Such that `feedback(sys, cont)` produces a closed-loop system with eigenvalues given by `A-KC` and `A-BL`.

# Arguments:
- `sys`: Model of system
- `L`: State-feedback gain `u = -Lx`
- `K`: Observer gain

See also [`observer_predictor`](@ref) and [`innovation_form`](@ref).
"""
function observer_controller(sys, L::AbstractMatrix, K::AbstractMatrix)
    A,B,C,D = ssdata(sys)
    ss(A - B*L - K*C + K*D*L, K, L, 0, sys.timeevol)
end
