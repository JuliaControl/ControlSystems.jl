"""`care(A, B, Q, R)`

Compute 'X', the solution to the continuous-time algebraic Riccati equation,
defined as A'X + XA - (XB)R^-1(B'X) + Q = 0, where R is non-singular.

Algorithm taken from:
Laub, "A Schur Method for Solving Algebraic Riccati Equations."
http://dspace.mit.edu/bitstream/handle/1721.1/1301/R-0859-05666488.pdf
"""
function care(A, B, Q, R)
    G = try
        B*inv(R)*B'
    catch y
        if y isa SingularException
            error("R must be non-singular in care.")
        else
            throw(y)
        end
    end

    Z = [A  -G;
        -Q  -A']

    S = schur(Z)
    S = ordschur(S, real(S.values).<0)
    U = S.Z

    (m, n) = size(U)
    U11 = U[1:div(m, 2), 1:div(n,2)]
    U21 = U[div(m,2)+1:m, 1:div(n,2)]
    return U21/U11
end

"""`dare(A, B, Q, R)`

Compute `X`, the solution to the discrete-time algebraic Riccati equation,
defined as A'XA - X - (A'XB)(B'XB + R)^-1(B'XA) + Q = 0, where Q>=0 and R>0

Algorithm taken from:
Laub, "A Schur Method for Solving Algebraic Riccati Equations."
http://dspace.mit.edu/bitstream/handle/1721.1/1301/R-0859-05666488.pdf
"""
function dare(A, B, Q, R)
    if (!ishermitian(Q) || minimum(eigvals(real(Q))) < 0)
        error("Q must be positive-semidefinite.");
    end
    if (!isposdef(R))
        error("R must be positive definite.");
    end

    n = size(A, 1);

    E = [
        Matrix{Float64}(I, n, n) B/R*B';
        zeros(size(A)) A'
    ];
    F = [
        A zeros(size(A));
        -Q Matrix{Float64}(I, n, n)
    ];

    QZ = schur(F, E);
    QZ = ordschur(QZ, abs.(QZ.alpha./QZ.beta) .< 1);

    return QZ.Z[(n+1):end, 1:n]/QZ.Z[1:n, 1:n];
end

"""`dlyap(A, Q)`

Compute the solution `X` to the discrete Lyapunov equation
`AXA' - X + Q = 0`.
"""
function dlyap(A::AbstractMatrix, Q)
    lhs = kron(A, conj(A))
    lhs = I - lhs
    x = lhs\reshape(Q, prod(size(Q)), 1)
    return reshape(x, size(Q))
end

"""`gram(sys, opt)`

Compute the grammian of system `sys`. If `opt` is `:c`, computes the
controllability grammian. If `opt` is `:o`, computes the observability
grammian."""
function gram(sys::AbstractStateSpace, opt::Symbol)
    if !isstable(sys)
        error("gram only valid for stable A")
    end
    func = iscontinuous(sys) ? lyap : dlyap
    if opt == :c
        # TODO probably remove type check in julia 0.7.0
        return func(sys.A, sys.B*sys.B')#::Array{numeric_type(sys),2} # lyap is type-unstable
    elseif opt == :o
        return func(Matrix(sys.A'), sys.C'*sys.C)#::Array{numeric_type(sys),2} # lyap is type-unstable
    else
        error("opt must be either :c for controllability grammian, or :o for
                observability grammian")
    end
end

"""`obsv(A, C)` or `obsv(sys)`

Compute the observability matrix for the system described by `(A, C)` or `sys`.

Note that checking for observability by computing the rank from `obsv` is
not the most numerically accurate way, a better method is checking if
`gram(sys, :o)` is positive definite."""
function obsv(A::AbstractMatrix, C::AbstractMatrix)
    T = promote_type(eltype(A), eltype(C))
    n = size(A, 1)
    ny = size(C, 1)
    if n != size(C, 2)
        throw(ArgumentError("C must have the same number of columns as A"))
    end
    res = fill(zero(T), n*ny, n)
    res[1:ny, :] = C
    for i=1:n-1
        res[(1 + i*ny):(1 + i)*ny, :] = res[((i - 1)*ny + 1):i*ny, :] * A
    end
    return res
end
obsv(sys::StateSpace) = obsv(sys.A, sys.C)

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
ctrb(sys::StateSpace) = ctrb(sys.A, sys.B)

"""`P = covar(sys, W)`

Calculate the stationary covariance `P = E[y(t)y(t)']` of an lti-model `sys`, driven by gaussian
white noise 'w' of covariance `E[w(t)w(τ)]=W*δ(t-τ)` where δ is the dirac delta.

The ouput is if Inf if the system is unstable. Passing white noise directly to
the output will result in infinite covariance in the corresponding outputs
(D*W*D' .!= 0) for contunuous systems."""
function covar(sys::AbstractStateSpace, W)
    (A, B, C, D) = (sys.A, sys.B, sys.C, sys.D)
    if !isa(W, UniformScaling) && (size(B,2) != size(W, 1) || size(W, 1) != size(W, 2))
        error("W must be a square matrix the same size as `sys.B` columns")
    end
    if !isstable(sys)
        return fill(Inf,(size(C,1),size(C,1)))
    end
    func = iscontinuous(sys) ? lyap : dlyap
    Q = try
        func(A, B*W*B')
    catch
        error("No solution to the Lyapunov equation was found in covar")
    end
    P = C*Q*C'
    if iscontinuous(sys)
        #Variance and covariance infinite for direct terms
        directNoise = D*W*D'
        for i in 1:size(C,1)
            if directNoise[i,i] != 0
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

"""
    covar(C,W)
If `C` is a matrix, return CWC'
"""
covar(C::Union{AbstractMatrix,UniformScaling}, R) = C*R*C'


# Note: the H∞ norm computation is probably not as accurate as with SLICOT,
# but this seems to be still reasonably ok as a first step
"""
`..  norm(sys, p=2; tol=1e-6)`

`norm(sys)` or `norm(sys,2)` computes the H2 norm of the LTI system `sys`.

`norm(sys, Inf)` computes the L∞ norm of the LTI system `sys`.
The H∞ norm is the same as the L∞ for stable systems, and Inf for unstable systems.
If the peak gain frequency is required as well, use the function `norminf` instead.

`tol` is an optional keyword argument, used only for the computation of L∞ norms.
It represents the desired relative accuracy for the computed L∞ norm
(this is not an absolute certificate however).

sys is first converted to a state space model if needed.

The L∞ norm computation implements the 'two-step algorithm' in:
N.A. Bruinsma and M. Steinbuch, 'A fast algorithm to compute the H∞-norm
of a transfer function matrix', Systems and Control Letters 14 (1990), pp. 287-293.
For the discrete-time version, see, e.g.,: P. Bongers, O. Bosgra, M. Steinbuch, 'L∞-norm
calculation for generalized state space systems in continuous and discrete time',
American Control Conference, 1991.
"""
function LinearAlgebra.norm(sys::AbstractStateSpace, p::Real=2; tol=1e-6)
    if p == 2
        return sqrt(tr(covar(sys, I)))
    elseif p == Inf
        if sys.Ts == 0
            return normLinf_twoSteps_ct(sys,tol)[1]
        else
            return normLinf_twoSteps_dt(sys,tol)[1]
        end
    else
        error("`p` must be either `2` or `Inf`")
    end
end

function LinearAlgebra.norm(sys::TransferFunction, p::Real=2; tol=1e-6)
    return norm(ss(sys), p, tol=tol)
end

"""
`.. (peakgain, peakgainfrequency) = norminf(sys; tol=1e-6)`

Compute the L∞ norm of the LTI system `sys`, together with the frequency
`peakgainfrequency` (in rad/TimeUnit) at which the gain achieves its peak value `peakgain`.
The H∞ norm is the same as the L∞ for stable systems, and Inf for unstable systems.

`tol` is an optional keyword argument representing the desired relative accuracy for
the computed L∞ norm (this is not an absolute certificate however).

sys is first converted to a state space model if needed.

The L∞ norm computation implements the 'two-step algorithm' in:
N.A. Bruinsma and M. Steinbuch, 'A fast algorithm to compute the H∞-norm
of a transfer function matrix', Systems and Control Letters 14 (1990), pp. 287-293.
For the discrete-time version, see, e.g.,: P. Bongers, O. Bosgra, M. Steinbuch, 'L∞-norm
calculation for generalized state space systems in continuous and discrete time',
American Control Conference, 1991.
"""
function norminf(sys::AbstractStateSpace; tol=1e-6)
    if sys.Ts == 0
        return normLinf_twoSteps_ct(sys,tol)
    else
        return normLinf_twoSteps_dt(sys,tol)
    end
end

function norminf(sys::TransferFunction, ; tol=1e-6)
    return norminf(ss(sys), tol=tol)
end

function normLinf_twoSteps_ct(sys::AbstractStateSpace, tol=1e-6, maxIters=1000, approximag=1e-10)
    # `maxIters`: the maximum  number of iterations allowed in the algorithm (default 1000)
    # approximag is a tuning parameter: what does it mean for a number to be on the imaginary axis
    # Because of this tuning for example, the relative precision that we provide on the norm computation
    # is not a true guarantee, more an order of magnitude
    # outputs: pair of Float64, namely L∞ norm approximation and frequency fpeak at which it is achieved
    T = promote_type(numeric_type(sys), Float64)
    if sys.nx == 0  # static gain
        return (svdvals(sys.D)[1], T(0))
    end
    p = pole(sys)
    # Check if there is a pole on the imaginary axis
    pidx = findfirst(map(x->isapprox(x,0.0),real(p)))
    if !(pidx isa Nothing)
        return (T(Inf), imag(p[pidx]))
        # note: in case of cancellation, for s/s for example, we return Inf, whereas Matlab returns 1
    else
        # Initialization: computation of a lower bound from 3 terms
        lb = maximum(svdvals(sys.D))
        fpeak = T(Inf)
        (lb, idx) = findmax([lb, T(maximum((svdvals(evalfr(sys,0)))))]) #TODO remove T() in julia 0.7.0
        if idx == 2
            fpeak = T(0)
        end
        if isreal(p)  # only real poles
            omegap = minimum(abs.(p))
        else  # at least one pair of complex poles
            tmp = abs.(imag.(p)./(real.(p).*abs.(p)))
            omegap = abs(p[argmax(tmp)])    # TODO This is highly suspicious
        end
        (lb, idx) = findmax([lb, T(maximum(svdvals(evalfr(sys, omegap*1im))))]) #TODO remove T() in julia 0.7.0
        if idx == 2
            fpeak = omegap
        end

        # Iterations
        iter = 1;
        while iter <= maxIters
            res = (1+2*T(tol))*lb
            R = sys.D'*sys.D - res^2*I
            S = sys.D*sys.D' - res^2*I
            M = sys.A-sys.B*(R\sys.D')*sys.C
            H = [         M              -res*sys.B*(R\sys.B') ;
                   res*sys.C'*(S\sys.C)            -M'            ]
            omegas = eigvals(H) .+ 0im # To make type stable
            omegaps = imag.(omegas[ (abs.(real.(omegas)).<=approximag) .& (imag.(omegas).>=0) ])
            sort!(omegaps)
            if isempty(omegaps)
                return (1+T(tol))*lb, fpeak
            else  # if not empty, omegaps contains at least two values
                ms = [(x+y)/2 for x=omegaps[1:end-1], y=omegaps[2:end]]
                for mval in ms
                    (lb, idx) = findmax([lb, T(maximum(svdvals(evalfr(sys,mval*1im))))]) #TODO remove T() in julia 0.7.0
                    if idx == 2
                        fpeak = mval
                    end
                end
            end
            iter += 1
        end
        error("In norminf: The computation of the H-infinity norm did not converge in $maxIters iterations")
    end
end

# discrete-time version of normHinf_twoSteps_ct above
# The value fpeak returned by the function is in the range [0,pi)/sys.Ts (in rad/s)
function normLinf_twoSteps_dt(sys::AbstractStateSpace,tol=1e-6, maxIters=1000, approxcirc=1e-8)
    T = promote_type(numeric_type(sys), Float64)
    if sys.nx == 0  # static gain
        return (svdvals(sys.D)[1], T(0))
    end
    p = pole(sys)
    # Check first if there is a pole on the unit circle
    pidx = findfirst(map(x->isapprox(x,1.0),abs.(p)))
    if !(pidx isa Nothing)
        return (T(Inf), angle(p[pidx])/abs(T(sys.Ts)))
    else
        # Initialization: computation of a lower bound from 3 terms
        lb = T(maximum(svdvals(evalfr(sys,1))))  #TODO remove T() in julia 0.7.0
        fpeak = T(0)
        (lb, idx) = findmax([lb, T(maximum(svdvals(evalfr(sys,-1))))]) #TODO remove T() in julia 0.7.0
        if idx == 2
            fpeak = T(pi)
        end

        p = p[imag(p).>0]
        if ~isempty(p)  # not just real poles
            # find frequency of pôle closest to unit circle
            omegap = angle(p[findmin(abs.(abs.(p).-1))[2]])
        else
            omegap = T(pi)/2
        end
        (lb, idx) = findmax([lb, T(maximum(svdvals(evalfr(sys, exp(omegap*1im)))))]) #TODO remove T() in julia 0.7.0
        if idx == 2
            fpeak = omegap
        end

        # Iterations
        iter = 1;
        while iter <= maxIters
            res = (1+2*T(tol))*lb
            R = res^2*I - sys.D'*sys.D
            RinvDt = R\sys.D'
            L = [ sys.A+sys.B*RinvDt*sys.C  sys.B*(R\sys.B');
                  zeros(T, sys.nx,sys.nx)      I]
            M = [ I                                 zeros(T, sys.nx,sys.nx);
                  sys.C'*(I+sys.D*RinvDt)*sys.C     L[1:sys.nx,1:sys.nx]']
            # +0im to make type stable
            zs = eigvals(L,M) .+ 0im # generalized eigenvalues
            # are there eigenvalues on the unit circle?
            omegaps = angle.(zs[ (abs.(abs.(zs).-1) .<= approxcirc) .& (imag(zs).>=0)])
            sort!(omegaps)
            if isempty(omegaps)
                return (1+T(tol))*lb, fpeak/T(sys.Ts)
            else  # if not empty, omegaps contains at least two values
                ms = [(x+y)/2 for x=omegaps[1:end-1], y=omegaps[2:end]]
                for mval in ms
                    (lb, idx) = findmax([lb, T(maximum(svdvals(evalfr(sys,exp(mval*1im)))))]) #TODO remove T() in julia 0.7.0
                    if idx == 2
                        fpeak = mval
                    end
                end
            end
            iter += 1
        end
        error("In norminf: The computation of the H-infinity norm did not converge in $maxIters iterations")
    end
end


"""`S, P, B = balance(A[, perm=true])`

Compute a similarity transform `T` resulting in `B = T\\A*T` such that the row
and column norms of `B` are approximately equivalent. If `perm=false`, the
transformation will only scale `A` using diagonal `S`, and not permute `A` (i.e., set `P=I`)."""
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
`sysr, G = balreal(sys::StateSpace)`

Calculates a balanced realization of the system sys, such that the observability and reachability gramians of the balanced system are equal and diagonal `G`

See also `gram`, `baltrunc`

Glad, Ljung, Reglerteori: Flervariabla och Olinjära metoder
"""
function balreal(sys::ST) where ST <: AbstractStateSpace
    P = gram(sys, :c)
    Q = gram(sys, :o)

    Q1 = try
        cholesky(Hermitian(Q)).U
    catch
        throw(ArgumentError("Balanced realization failed: Observability grammian not positive definite, system needs to be observable"))
    end
    U,Σ,V = svd(Q1*P*Q1')
    Σ .= sqrt.(Σ)
    Σ1 = diagm(0 => sqrt.(Σ))
    T = Σ1\(U'Q1)

    Pz = T*P*T'
    Qz = inv(T')*Q*inv(T)
    if norm(Pz-Qz) > sqrt(eps())
        @warn("balreal: Result may be inaccurate")
        println("Controllability gramian before transform")
        display(P)
        println("Controllability gramian after transform")
        display(Pz)
        println("Observability gramian before transform")
        display(Q)
        println("Observability gramian after transform")
        display(Qz)
        println("Singular values of PQ")
        display(Σ)
    end

    sysr = ST(T*sys.A/T, T*sys.B, sys.C/T, sys.D, sys.Ts), diagm(0 => Σ)
end


"""
`sysr, G = baltrunc(sys::StateSpace, atol = √ϵ, rtol=1e-3, unitgain=true)`

Reduces the state dimension by calculating a balanced realization of the system sys, such that the observability and reachability gramians of the balanced system are equal and diagonal `G`, and truncating it such that all states corresponding to singular values less than `atol` and less that `rtol σmax` are removed. If `unitgain=true`, the matrix `D` is chosen such that unit static gain is achieved.

See also `gram`, `balreal`

Glad, Ljung, Reglerteori: Flervariabla och Olinjära metoder
"""
function baltrunc(sys::ST; atol = sqrt(eps()), rtol = 1e-3, unitgain = true) where ST <: AbstractStateSpace
    sysbal, S = balreal(sys)
    S = diag(S)
    S = S[S .>= atol]
    S = S[S .>= S[1]*rtol]
    n = length(S)
    A = sysbal.A[1:n,1:n]
    B = sysbal.B[1:n,:]
    C = sysbal.C[:,1:n]
    D = sysbal.D
    if unitgain
        D = D/(C*inv(-A)*B)
    end

    return ST(A,B,C,D,sys.Ts), diagm(0 => S)
end

"""
    syst = similarity_transform(sys, T)
Perform a similarity transform `T : Tx̃ = x` on `sys` such that
```
Ã = T⁻¹AT
B̃ = T⁻¹ B
C̃ = CT
D̃ = D
```
"""
function similarity_transform(sys::ST, T) where ST <: AbstractStateSpace
    Tf = factorize(T)
    A = Tf\sys.A*T
    B = Tf\sys.B
    C = sys.C*T
    D = sys.D
    ST(A,B,C,D,sys.Ts)
end

"""
sysi = innovation_form(sys, R1, R2)
sysi = innovation_form(sys; sysw=I, syse=I, R1=I, R2=I)

Takes a system
```
x' = Ax + Bu + w ~ R1
y  = Cx + e ~ R2
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
function innovation_form(sys::ST, R1, R2) where ST <: AbstractStateSpace
    K = kalman(sys, R1, R2)
    ST(sys.A, K, sys.C, Matrix{eltype(sys.A)}(I, sys.ny, sys.ny), sys.Ts)
end
# Set D = I to get transfer function H = I + C(sI-A)\ K
function innovation_form(sys::ST; sysw=I, syse=I, R1=I, R2=I) where ST <: AbstractStateSpace
	K = kalman(sys, covar(sysw,R1), covar(syse, R2))
	ST(sys.A, K, sys.C, Matrix{eltype(sys.A)}(I, sys.ny, sys.ny), sys.Ts)
end
