import Base.LinAlg: BlasFloat

@doc """`care(A, B, Q, R)`

Compute 'X', the solution to the continuous-time algebraic Riccati equation,
defined as A'X + XA - (XB)R^-1(B'X) + Q = 0, where R is non-singular.

Algorithm taken from:
Laub, "A Schur Method for Solving Algebraic Riccati Equations."
http://dspace.mit.edu/bitstream/handle/1721.1/1301/R-0859-05666488.pdf
""" ->
function care(A, B, Q, R)
    G = try
        B*inv(R)*B'
    catch
        error("R must be non-singular.")
    end

    Z = [A  -G;
        -Q  -A']

    S = schurfact(Z)
    S = ordschur(S, real(S.values).<0)
    U = S.Z

    (m, n) = size(U)
    U11 = U[1:div(m, 2), 1:div(n,2)]
    U21 = U[div(m,2)+1:m, 1:div(n,2)]
    return U21/U11
end

@doc """`dare(A, B, Q, R)`

Compute `X`, the solution to the discrete-time algebraic Riccati equation,
defined as A'XA - X - (A'XB)(B'XB + R)^-1(B'XA) + Q = 0, where A and R
are non-singular.

Algorithm taken from:
Laub, "A Schur Method for Solving Algebraic Riccati Equations."
http://dspace.mit.edu/bitstream/handle/1721.1/1301/R-0859-05666488.pdf
""" ->
function dare(A, B, Q, R)
    G = try
        B*inv(R)*B'
    catch
        error("R must be non-singular.")
    end

    Ait = try
        inv(A)'
    catch
        error("A must be non-singular.")
    end

    Z = [A + G*Ait*Q   -G*Ait;
         -Ait*Q        Ait]

    S = schurfact(Z)
    S = ordschur(S, abs(S.values).<=1)
    U = S.Z

    (m, n) = size(U)
    U11 = U[1:div(m, 2), 1:div(n,2)]
    U21 = U[div(m,2)+1:m, 1:div(n,2)]
    return U21/U11
end

@doc """`dlyap(A, Q)`

Compute the solution "X" to the discrete Lyapunov equation
"AXA' - X + Q = 0".
"""
function dlyap(A, Q)
    lhs = kron(A, conj(A))
    lhs = eye(size(lhs, 1)) - lhs
    x = lhs\reshape(Q, prod(size(Q)), 1)
    return reshape(x, size(Q))
end

@doc """`gram(sys, opt)`

Compute the grammian of system `sys`. If `opt` is `:c`, computes the
controllability grammian. If `opt` is `:o`, computes the observability
grammian.""" ->
function gram(sys::StateSpace, opt::Symbol)
    if !isstable(sys)
        error("gram only valid for stable A")
    end
    func = iscontinuous(sys) ? lyap : dlyap
    if opt == :c
        return func(sys.A, sys.B*sys.B')
    elseif opt == :o
        return func(sys.A', sys.C'*sys.C)
    else
        error("opt must be either :c for controllability grammian, or :o for
                observability grammian")
    end
end

@doc """`obsv(A, C)` or `obsv(sys)`

Compute the observability matrix for the system described by `(A, C)` or `sys`.

Note that checking for observability by computing the rank from `obsv` is
not the most numerically accurate way, a better method is checking if
`gram(sys, :o)` is positive definite.""" ->
function obsv(A, C)
    n = size(A, 1)
    ny = size(C, 1)
    if n != size(C, 2)
        error("C must have the same number of columns as A")
    end
    res = zeros(n*ny, n)
    res[1:ny, :] = C
    for i=1:n-1
        res[(1 + i*ny):(1 + i)*ny, :] = res[((i - 1)*ny + 1):i*ny, :] * A
    end
    return res
end
obsv(sys::StateSpace) = obsv(sys.A, sys.C)

@doc """`ctrb(A, B)` or `ctrb(sys)`

Compute the controllability matrix for the system described by `(A, B)` or
`sys`. 

Note that checking for controllability by computing the rank from
`obsv` is not the most numerically accurate way, a better method is
checking if `gram(sys, :c)` is positive definite.""" ->
function ctrb(A, B)
    n = size(A, 1)
    nu = size(B, 2)
    if n != size(B, 1)
        error("B must have the same number of rows as A")
    end
    res = zeros(n, n*nu)
    res[:, 1:nu] = B
    for i=1:n-1
        res[:, (1 + i*nu):(1 + i)*nu] = A * res[:, ((i - 1)*nu + 1):i*nu]
    end
    return res
end
ctrb(sys::StateSpace) = ctrb(sys.A, sys.B)

@doc """`covar(sys, W)`

Calculate the stationary covariance of an lti-model `sys`, driven by gaussian
white noise of covariance `W`""" ->
function covar(sys::StateSpace, W::StridedMatrix)
    (A, B, C, D) = (sys.A, sys.B, sys.C, sys.D)
    if size(B,2) != size(W, 1) || size(W, 1) != size(W, 2)
        error("W must be a square matrix the same size as `sys.B` columns")
    end
    if !isstable(sys) || any(D .!= 0)
        return Inf
    end
    func = iscontinuous(sys) ? lyap : dlyap
    Q = func(A, B*W*B')
    return C*Q*C' + D*W*D'
end

@doc """`norm(sys[, p])`

Compute the `p`-norm of the system `sys`. `p` can be either `2` or `Inf`
(default is 2)""" ->
function Base.norm(sys::StateSpace, p::Real=2)
    if p == 2
        return sqrt(trace(covar(sys, eye(size(sys.B, 2)))))
    elseif p == Inf
        error("Hinf norm not implemented")
    else
        error("`p` must be either `2` or `Inf`")
    end
end

@doc """`T, B = balance(A[, perm=true])`

Compute a similarity transform `T` resulting in `B = T\\A*T` such that the row
and column norms of `B` are approximately equivalent. If `perm=false`, the
transformation will only scale, and not permute `A`.""" ->
function balance(A, perm::Bool=true)
    n = Base.LinAlg.chksquare(A)
    B = copy(A)
    job = perm ? 'B' : 'S'
    ilo, ihi, scale = LAPACK.gebal!(job, B)

    S = diagm(scale)
    for j = 1:(ilo-1)   S[j,j] = 1 end
    for j = (ihi+1):n   S[j,j] = 1 end

    P = eye(Int, n)
    if perm
        if ilo > 1
            for j = (ilo-1):-1:1 cswap!(j, round(Int, scale[j]), P) end
        end
        if ihi < n
            for j = (ihi+1):n    cswap!(j, round(Int, scale[j]), P) end
        end
    end
    return S, P, B
end

function cswap!{T<:Number}(i::Integer, j::Integer, X::StridedMatrix{T})
    for k = 1:size(X,1)
        X[i, k], X[j, k] = X[j, k], X[i, k]
    end
end
