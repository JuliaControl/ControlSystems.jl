const _scaling_notice = """
Note: Gramian computations are sensitive to input-output scaling. For the result of a numerical balancing, gramian computation or truncation of MIMO systems to be meaningful, the inputs and outputs of the system must thus be scaled in a meaningful way. A common (but not the only) approach is:
- The outputs are scaled such that the maximum allowed control error, the maximum expected reference variation, or the maximum expected variation, is unity.
- The input variables are scaled to have magnitude one. This is done by dividing each variable by its maximum expected or allowed change, i.e., ``u_{scaled} = u / u_{max}``

Without such scaling, the result of balancing will depend on the units used to measure the input and output signals, e.g., a change of unit for one output from meter to millimeter will make this output 1000x more important.
"""

"""
    are(::Continuous, A, B, Q, R)

Compute 'X', the solution to the continuous-time algebraic Riccati equation,
defined as A'X + XA - (XB)R^-1(B'X) + Q = 0, where R is non-singular.

In an LQR problem, `Q` is associated with the state penalty ``x'Qx`` while `R` is associated with the control penalty ``u'Ru``.
See [`lqr`](@ref) for more details.

Uses `MatrixEquations.arec`. For keyword arguments, see the docstring of `ControlSystemsBase.MatrixEquations.arec`,
note that they define the input arguments in a different order.
"""
function are(::ContinuousType, A::AbstractMatrix, B, Q, R, args...; kwargs...)
    arec(A, B, R, Q, args...; kwargs...)[1]
end

"""
    are(::Discrete, A, B, Q, R; kwargs...)

Compute `X`, the solution to the discrete-time algebraic Riccati equation,
defined as A'XA - X - (A'XB)(B'XB + R)^-1(B'XA) + Q = 0, where Q>=0 and R>0

In an LQR problem, `Q` is associated with the state penalty ``x'Qx`` while `R` is associated with the control penalty ``u'Ru``.
See [`lqr`](@ref) for more details.

Uses `MatrixEquations.ared`. For keyword arguments, see the docstring of `ControlSystemsBase.MatrixEquations.ared`,
note that they define the input arguments in a different order.
"""
function are(::DiscreteType, A::AbstractMatrix, B, Q, R, args...; kwargs...)
    ared(A, B, R, Q, args...; kwargs...)[1]
end

are(t::TimeEvolType, A::Number, B::Number, Q::Number, R::Number) = are(t, fill(A,1,1),fill(B,1,1),fill(Q,1,1),fill(R,1,1))
are(sys::AbstractStateSpace, args...; kwargs...) = are(timeevol(sys), sys.A, sys.B, args...; kwargs...)

@deprecate care(args...; kwargs...) are(Continuous, args...; kwargs...)
@deprecate dare(args...; kwargs...) are(Discrete, args...; kwargs...)

"""
    lyap(A, Q; kwargs...)

Compute the solution `X` to the discrete Lyapunov equation
`AXA' - X + Q = 0`.

Uses `MatrixEquations.lyapc / MatrixEquations.lyapd`. For keyword arguments, see the docstring of `ControlSystemsBase.MatrixEquations.lyapc / ControlSystemsBase.MatrixEquations.lyapd`
"""
function LinearAlgebra.lyap(::DiscreteType, A::AbstractMatrix, Q; kwargs...)
    lyapd(A, Q; kwargs...)
end

LinearAlgebra.lyap(::ContinuousType, args...; kwargs...) = lyapc(args...; kwargs...)
LinearAlgebra.lyap(::DiscreteType, args...; kwargs...) = lyapd(args...; kwargs...)
LinearAlgebra.lyap(sys::AbstractStateSpace, args...; kwargs...) = lyap(timeevol(sys), sys.A, args...; kwargs...)

"""
    Xc = plyap(sys::AbstractStateSpace, Ql; kwargs...)

Lyapunov solver that takes the `L` Cholesky factor of `Q` and returns a triangular matrix `Xc` such that `Xc*Xc' = X`.
"""
plyap(sys::AbstractStateSpace, args...; kwargs...) = plyap(timeevol(sys), sys.A, args...; kwargs...)
plyap(::ContinuousType, args...; kwargs...) = MatrixEquations.plyapc(args...; kwargs...)
plyap(::DiscreteType, args...; kwargs...) = MatrixEquations.plyapd(args...; kwargs...)

@deprecate dlyap(args...; kwargs...) lyap(Discrete, args...; kwargs...)


"""
    U = grampd(sys, opt; kwargs...)

Return a Cholesky factor `U` of the grammian of system `sys`. If `opt` is `:c`, computes the
controllability grammian `G = U*U'`. If `opt` is `:o`, computes the observability
grammian `G = U'U`.

Obtain a `Cholesky` object by `Cholesky(U)` for observability grammian

Uses `MatrixEquations.plyapc/plyapd`. For keyword arguments, see the docstring of `ControlSystemsBase.MatrixEquations.plyapc/plyapd`
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

# Extended help
$(_scaling_notice)
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
`gram(sys, :o)` is positive definite or to call the function [`observability`](@ref).

The unobservable subspace is `nullspace(obsv(A, C))`, initial conditions in this subspace produce a zero response.
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

"""
    ctrb(A, B)
    ctrb(sys)

Compute the controllability matrix for the system described by `(A, B)` or
`sys`.

Note that checking for controllability by computing the rank from
`ctrb` is not the most numerically accurate way, a better method is
checking if `gram(sys, :c)` is positive definite or to call the function [`controllability`](@ref).

The controllable subspace is given by the range of this matrix, and the uncontrollable subspace is `nullspace(ctrb(A, B)') (note the transpose)`.
"""
function ctrb(A::AbstractMatrix, B::AbstractVecOrMat)
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

"""
    controllability(A, B; atol, rtol)
    controllability(sys; atol, rtol)

Check for controllability of the pair `(A, B)` or `sys` using the PHB test.

The return value contains the field `iscontrollable` which is `true` if the rank condition is met at all eigenvalues of `A`, and `false` otherwise. The returned structure also contains the rank and smallest singular value at each individual eigenvalue of `A` in the fields `ranks` and `sigma_min`.

Technically, this function checks for controllability from the origin, also called reachability.
"""
function controllability(A::AbstractMatrix{T}, B; atol::Real=0, rtol::Real=atol>0 ? 0 : size(A,1)*eps(float(T))) where T
    n = LinearAlgebra.checksquare(A)
    p = eigvals(A)
    ranks = zeros(Int, n)
    sigma_min = similar(A, float(T), n)
    for i = 1:n
        sigmas = svdvals([(p[i]*I - A) B])
        r = count(>=(max(atol, rtol*sigmas[1])), sigmas)
        ranks[i] = r
        sigma_min[i] = sigmas[end]
    end
    (; iscontrollable = all(==(n), ranks), ranks, sigma_min)
end
controllability(sys::AbstractStateSpace; kwargs...) = controllability(sys.A, sys.B; kwargs...)


"""
    observability(A, C; atol, rtol)


Check for observability of the pair `(A, C)` or `sys` using the PHB test.

The return value contains the field `isobservable` which is `true` if the rank condition is met at all eigenvalues of `A`, and `false` otherwise. The returned structure also contains the rank and smallest singular value at each individual eigenvalue of `A` in the fields `ranks` and `sigma_min`.
"""
function observability(A::AbstractMatrix{T}, C; atol::Real=0, rtol::Real=atol>0 ? 0 : size(A,1)*eps(float(T))) where T
    n = LinearAlgebra.checksquare(A)
    p = eigvals(A)
    ranks = zeros(Int, n)
    sigma_min = similar(A, float(T), n)
    for i = 1:n
        sigmas = svdvals([(p[i]*I - A); C])
        r = count(>=(max(atol, rtol*sigmas[1])), sigmas)
        ranks[i] = r
        sigma_min[i] = sigmas[end]
    end
    (; isobservable = all(==(n), ranks), ranks, sigma_min)
end
observability(sys::AbstractStateSpace; kwargs...) = observability(sys.A, sys.C; kwargs...)

"""
    P = covar(sys, W)

Calculate the stationary covariance `P = E[y(t)y(t)']` of the output `y` of a
`StateSpace` model `sys` driven by white Gaussian noise `w` with covariance
`E[w(t)w(τ)]=W*δ(t-τ)` (δ is the Dirac delta).

Remark: If `sys` is unstable then the resulting covariance is a matrix of `Inf`s.
Entries corresponding to direct feedthrough (D*W*D' .!= 0) will equal `Inf`
for continuous-time systems.
    
See also [`innovation_form`](@ref).
"""
function covar(sys::AbstractStateSpace, W)
    (A, B, C, D) = ssdata(sys)
    if !isa(W, UniformScaling) && (size(B,2) != size(W, 1) || size(W, 1) != size(W, 2))
        error("W must be a square matrix the same size as `sys.B` columns")
    end
    isa(W, UniformScaling) && (W = W(size(B, 2)))
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
The H∞ norm is the same as the L∞ for stable systems, and Inf for unstable systems.
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
    sysm, T, SF = schur_form(sys)

Bring `sys` to Schur form.

The Schur form is characterized by `A` being Schur with the real values of eigenvalues of `A` on the main diagonal. `T` is the similarity transform applied to the system such that 
```julia
sysm ≈ similarity_transform(sys, T)
```
`SF` is the Schur-factorization of `A`.
"""
function schur_form(sys)
    SF = schur(sys.A)
    A = SF.T
    B = SF.Z'*sys.B
    C = sys.C*SF.Z
    ss(A,B,C,sys.D, sys.timeevol), SF.Z, SF
end

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
hinfnorm(sys::AbstractStateSpace{<:Continuous}; tol=1e-6) = _infnorm_two_steps_ct(schur_form(sys)[1], :hinf, tol)
hinfnorm(sys::AbstractStateSpace{<:Discrete}; tol=1e-6) = _infnorm_two_steps_dt(schur_form(sys)[1], :hinf, tol)
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
    sys2, _ = schur_form(sys)
    if iscontinuous(sys2)
        return _infnorm_two_steps_ct(sys2, :linf, tol)
    else
        return _infnorm_two_steps_dt(sys2, :linf, tol)
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
    # outputs: An approximation of the L∞ norm and the frequency ω_peak at which it is achieved
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
    @error("In _infnorm_two_steps_dt: The computation of the H∞/L∞ norm did not converge in $maxIters iterations")
    return T((1+tol)*lb), T(ω_peak)
end

function _infnorm_two_steps_dt(sys::AbstractStateSpace, normtype::Symbol, tol=1e-6, maxIters=250, approxcirc=1e-8)
    # Discrete-time version of linfnorm_two_steps_ct above
    # Computations are done in normalized frequency θ

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
function balance(A::AbstractMatrix{<:LinearAlgebra.BlasFloat}, perm::Bool=true)
    n = LinearAlgebra.checksquare(A)
    B = copy(A)
    job = perm ? 'B' : 'S'
    ilo, ihi, scaling = LAPACK.gebal!(job, B)

    S = Diagonal(scaling)
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

function balance(A::AbstractMatrix, perm::Bool=true)
    Ac = Float64.(A)
    balance(Ac, perm)
end



"""
`sysr, G, T = balreal(sys::StateSpace)`

Calculates a balanced realization of the system sys, such that the observability and reachability gramians of the balanced system are equal and diagonal `diagm(G)`. `T` is the similarity transform between the old state `x` and the new state `z` such that `z = Tx`.

See also [`gram`](@ref), [`baltrunc`](@ref).

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
    rmin = count(Σ .> sqrt(eps(eltype(Σ)))*Σ[1])
    i1 = 1:rmin
    Σ = Σ[i1]


    hsi2 = Diagonal(1 ./sqrt.(Σ))
    L = lmul!(R',view(U,:,i1))*hsi2
    Tr = lmul!(S,V[:,i1])*hsi2
    # return the minimal balanced system
    T = L'SF.Z'
    return ss(L'SF.T*Tr, L'bs, cs*Tr, sys.D, sys.timeevol), Σ, T
end


"""
    sysr, G, T = baltrunc(sys::StateSpace; atol = √ϵ, rtol=1e-3, n = nothing, residual = false)

Reduces the state dimension by calculating a balanced realization of the system sys, such that the observability and reachability gramians of the balanced system are equal and diagonal `diagm(G)`, and truncating it to order `n`. If `n` is not provided, it's chosen such that all states corresponding to singular values less than `atol` and less that `rtol σmax` are removed.

`T` is the projection matrix between the old state `x` and the newstate `z` such that `z = Tx`. `T` will in general be a non-square matrix.

If `residual = true`, matched static gain is achieved through "residualization", i.e., setting
```math
0 = A_{21}x_{1} + A_{22}x_{2} + B_{2}u
```
where indices 1/2 correspond to the remaining/truncated states respectively.

See also `gram`, `balreal`

Glad, Ljung, Reglerteori: Flervariabla och Olinjära metoder.

For more advanced model reduction, see [RobustAndOptimalControl.jl - Model Reduction](https://juliacontrol.github.io/RobustAndOptimalControl.jl/dev/#Model-reduction).

# Extended help
$(_scaling_notice)
"""
function baltrunc(sys::ST; atol = sqrt(eps(numeric_type(sys))), rtol = 1e-3, n = nothing, residual=false) where ST <: AbstractStateSpace
    sysbal, S, T = balreal(sys)
    if n === nothing
        S = S[S .>= atol]
        S = S[S .>= S[1]*rtol]
        n = length(S)
    else
        n > sys.nx && error("n too large. A state dimension of n = $n was requested, but the original system has a $(sys.nx)-dimensional state.")
        if length(S) < n
            @error("n too large. A state dimension of n = $n was requested, but after a balanced realization was computed only $(length(S)) dimensions remain. Try either calling `minreal` before calling `baltrunc`, or try balancing the model using `balance_statespace`. Returning a system with n = $(length(S))")
            n = length(S)
        end
        S = S[1:n]
    end
    i1 = 1:n
    T = T[i1, :]
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

    return ST(A,B,C,D,sys.timeevol), S, T
end

"""
    stab, unstab, sep = stab_unstab(sys; kwargs...)

Decompose `sys` into `sys = stab + unstab` where `stab` contains all stable poles and `unstab` contains unstable poles. 

`0 ≤ sep ≤ 1` is the estimated separation between the stable and unstable spectra.

The docstring of `MatrixPencils.ssblkdiag`, reproduced below, provides more information on the keyword arguments:
$(@doc(MatrixPencils.ssblkdiag))
"""
function stab_unstab(sys::AbstractStateSpace; kwargs...)
    stable_unstable = true
    disc = isdiscrete(sys)
    A, B, C, _, _, blkdims, sep = MatrixPencils.ssblkdiag(sys.A, sys.B, sys.C; disc, stable_unstable, withQ = false, withZ = false, kwargs...)
    n1 = blkdims[1];
    i1 = 1:n1; i2 = n1+1:sys.nx 
    return (; stab=ss(A[i1,i1], B[i1,:], C[:,i1], sys.D, timeevol(sys)), 
            unstab=ss(A[i2,i2], B[i2,:], C[:,i2], 0, timeevol(sys)),
            sep)
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
See also [`balance_statespace`](@ref), [`find_similarity_transform`](@ref).
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
    find_similarity_transform(sys1, sys2, method = :obsv)

Find T such that `similarity_transform(sys1, T) == sys2`

Ref: Minimal state-space realization in linear system theory: an overview, B. De Schutter

If `method == :obsv`, the observability matrices of `sys1` and `sys2` are used to find `T`, whereas `method == :ctrb` uses the controllability matrices.

```jldoctest
julia> using ControlSystemsBase

julia> T = randn(3,3);

julia> sys1 = ssrand(1,1,3);

julia> sys2 = similarity_transform(sys1, T);

julia> T2 = find_similarity_transform(sys1, sys2);

julia> T2 ≈ T
true
```
"""
function find_similarity_transform(sys1, sys2, method = :obsv)
    if method === :obsv
        O1 = obsv(sys1)
        O2 = obsv(sys2)
        return O1\O2
    elseif method === :ctrb
        C1 = ctrb(sys1)
        C2 = ctrb(sys2)
        return C1/C2
    else
        error("Unknown method $method")
    end
end

"""
    time_scale(sys::AbstractStateSpace{Continuous}, a; balanced = false)
    time_scale(G::TransferFunction{Continuous},     a; balanced = true)

Rescale the time axis (change time unit) of `sys`.

For systems where the dominant time constants are very far from 1, e.g., in electronics, rescaling the time axis may be beneficial for numerical performance, in particular for continuous-time simulations.

Scaling of time for a function ``f(t)`` with Laplace transform ``F(s)`` can be stated as
```math
f(at) \\leftrightarrow \\dfrac{1}{a} F\\big(\\dfrac{s}{a}\\big)
```

The keyword argument `balanced` indicates whether or not to apply a balanced scaling on the `B` and `C` matrices.
For statespace systems, this defaults to false since it changes the state representation, only `B` will be scaled.
For transfer functions, it defaults to true.

# Example:
The following example show how a system with a time constant on the order of one micro-second is rescaled such that the time constant becomes 1, i.e., the time unit is changed from seconds to micro-seconds. 
```julia
Gs  = tf(1, [1e-6, 1])     # micro-second time scale modeled in seconds
Gms = time_scale(Gs, 1e-6) # Change to micro-second time scale
Gms == tf(1, [1, 1])       # Gms now has micro-seconds as time unit
```

The next example illustrates how the time axis of a time-domain simulation changes by time scaling 
```julia
t = 0:0.1:50 # original time axis
a = 10       # Scaling factor
sys1 = ssrand(1,1,5)
res1 = step(sys1, t)      # Perform original simulation
sys2 = time_scale(sys, a) # Scale time
res2 = step(sys2, t ./ a) # Simulate on scaled time axis, note the `1/a`
isapprox(res1.y, res2.y, rtol=1e-3, atol=1e-3)
```
"""
function time_scale(sys::AbstractStateSpace{Continuous}, a; balanced = false)
    a isa Real && (a > 0 || error("Time scaling constant must be positive"))
    A,B,C,D = ssdata(sys)
    A = a*I * A
    if balanced
        B = √(a)*I * B # Split the scaling equally on B and C to keep good balance
        C = √(a)*I * C
    else
        B = a*I*B
    end
    ss(A,B,C,D,sys.timeevol)
end

time_scale(G::TransferFunction{Continuous}, a; balanced = true) = tf(time_scale(ss(G), a; balanced))

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
	ss(sys.A, K, sys.C, Matrix{eltype(sys.A)}(I, sys.ny, sys.ny), sys.timeevol)
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
    observer_predictor(sys::AbstractStateSpace, K; h::Int = 1, output_state = false)
    observer_predictor(sys::AbstractStateSpace, R1, R2[, R12]; output_state = false)

If `sys` is continuous, return the observer predictor system
```math
\\begin{aligned}
x̂' &= (A - KC)x̂ + (B-KD)u + Ky \\\\
ŷ  &= Cx + Du
\\end{aligned}
```
with the input equation `[B-KD K] * [u; y]`

If `sys` is discrete, the prediction horizon `h` may be specified, in which case measurements up to and including time `t-h` and inputs up to and including time `t` are used to predict `y(t)`.

If covariance matrices `R1, R2` are given, the kalman gain `K` is calculated using [`kalman`](@ref).

If `output_state` is true, the output is the state estimate `x̂` instead of the output estimate `ŷ`.

See also [`innovation_form`](@ref), [`observer_controller`](@ref) and [`observer_filter`](@ref).
"""
function observer_predictor(sys::AbstractStateSpace, R1, R2::Union{AbstractArray, UniformScaling}, args...; kwargs...)
    K = kalman(sys, R1, R2, args...)
    observer_predictor(sys, K; kwargs...)
end

function observer_predictor(sys::AbstractStateSpace, K::AbstractMatrix; h::Integer = 1, output_state = false)
    h >= 1 || throw(ArgumentError("h must be positive."))
    ny = noutputs(sys)
    size(K, 1) == sys.nx && size(K,2) == ny || throw(ArgumentError("K has the wrong size, expected $((sys.nx, ny))"))
    A,B,C,D = ssdata(sys)
    if h == 1
        ss(A-K*C, [B-K*D K], output_state ? I : C, output_state ? 0 : [D zeros(ny, ny)], sys.timeevol)
    else
        isdiscrete(sys) || throw(ArgumentError("A prediction horizon is only supported for discrete systems. "))
        output_state && throw(ArgumentError("output_state is not supported for h > 1."))
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
    cont = observer_controller(sys, L::AbstractMatrix, K::AbstractMatrix; direct=false)


# If `direct = false`
Return the observer_controller `cont` that is given by
`ss(A - B*L - K*C + K*D*L, K, L, 0)`
such that `feedback(sys, cont)` produces a closed-loop system with eigenvalues given by `A-KC` and `A-BL`.

This controller does not have a direct term, and corresponds to state feedback operating on state estimated by [`observer_predictor`](@ref). Use this form if the computed control signal is applied at the next sampling instant, or with an otherwise large delay in relation to the measurement fed into the controller.

Ref: "Computer-Controlled Systems" Eq 4.37

# If `direct = true`
Return the observer controller `cont` that is given by
`ss((I-KC)(A-BL), (I-KC)(A-BL)K, L, LK)`
such that `feedback(sys, cont)` produces a closed-loop system with eigenvalues given by `A-BL` and `A-BL-KC`.
This controller has a direct term, and corresponds to state feedback operating on state estimated by [`observer_filter`](@ref). Use this form if the computed control signal is applied immediately after receiveing a measurement. This version typically has better performance than the one without a direct term.

!!! note
    To use this formulation, the observer gain `K` should have been designed for the pair `(A, CA)` rather than `(A, C)`. To do this, pass `direct = true` when calling [`place`](@ref) or [`kalman`](@ref).

Ref: Ref: "Computer-Controlled Systems" pp 140 and "Computer-Controlled Systems" pp 162 prob 4.7

# Arguments:
- `sys`: Model of system
- `L`: State-feedback gain `u = -Lx`
- `K`: Observer gain

See also [`observer_predictor`](@ref) and [`innovation_form`](@ref).
"""
function observer_controller(sys, L::AbstractMatrix, K::AbstractMatrix; direct=false)
    A,B,C,D = ssdata(sys)
    if direct && isdiscrete(sys)
        iszero(D) || throw(ArgumentError("D must be zero when using direct formulation of `observer_controller`"))
        IKC = (I - K*C)
        ABL = (A - B*L)
        ss(IKC*ABL,       IKC*ABL*K, L, L*K, sys.timeevol)
    else
        ss(A - B*L - K*C + K*D*L, K, L, 0,   sys.timeevol)
    end
end


"""
    observer_filter(sys, K; output_state = false)

Return the observer filter 
```math
\\begin{aligned}
x̂(k|k) &= (I - KC)Ax̂(k-1|k-1) + (I - KC)Bu(k-1) + Ky(k) \\\\
\\end{aligned}
```
with the input equation `[(I - KC)B K] * [u(k-1); y(k)]`.

Note the time indices in the equations, the filter assumes that the user passes the *current* ``y(k)``, but the *past* ``u(k-1)``, that is, this filter is used to estimate the state *before* the current control input has been applied. This causes a state-feedback controller acting on the estimate produced by this observer to have a direct term.

This is similar to [`observer_predictor`](@ref), but in contrast to the predictor, the filter output depends on the current measurement, whereas the predictor output only depend on past measurements.

The observer filter is equivalent to the [`observer_predictor`](@ref) for continuous-time systems.

!!! note
    To use this formulation, the observer gain `K` should have been designed for the pair `(A, CA)` rather than `(A, C)`. To do this, pass `direct = true` when calling [`place`](@ref) or [`kalman`](@ref).

Ref: "Computer-Controlled Systems" Eq 4.32
"""
function observer_filter(sys::AbstractStateSpace{<:Discrete}, K::AbstractMatrix; output_state = false)
    A,B,C,D = ssdata(sys)
    iszero(D) || throw(ArgumentError("D must be zero in `observer_filter`, consider using `observer_predictor` if you have a non-zero `D`."))
    IKC = (I-K*C)
    ss(IKC*A, [IKC*B K], output_state ? I : C, 0, sys.timeevol)
end

function observer_filter(sys::AbstractStateSpace{Continuous}, K::AbstractMatrix; kwargs...)
    observer_predictor(sys, K; kwargs...)
end