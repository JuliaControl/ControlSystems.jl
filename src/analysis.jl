@doc """`pole(sys)`

Compute the poles of system `sys`.""" ->
pole(sys::StateSpace) = eig(sys.A)[1]
pole(sys::TransferFunction) = [map(pole, sys.matrix)...]
pole(sys::SisoTf) = roots(sys.den)

@doc """`Wn, zeta = damp(sys)`

Compute the natural frequencies, `Wn`, and damping ratios, `zeta`, of the poles
of `sys`""" ->
function damp(sys::LTISystem)
    ps = pole(sys)
    if !iscontinuous(sys)
        Ts = sys.Ts == -1 ? 1 : sys.Ts
        ps = log(ps)/Ts
    end
    Wn = abs(ps)
    zeta = -cos(angle(ps))
    return Wn, zeta
end

@doc """`dampreport(sys)`

Display a report of the poles, damping ratio, natural frequency, and time
constant of the system `sys`""" ->
function dampreport(io::IO, sys::LTISystem)
    ps = pole(sys)
    np = length(ps)
    Wn, zeta = damp(sys)
    t_const = 1./(Wn.*zeta)
    order = sortperm(Wn)
    header =
    ("|     Pole      |   Damping     |   Frequency   | Time Constant |\n"*
     "|               |    Ratio      |   (rad/sec)   |     (sec)     |\n"*
     "+---------------+---------------+---------------+---------------+")
    println(io, header)
    for i=order
        p, z, w, t = ps[i], zeta[i], Wn[i], t_const[i]
        @printf(io, "|  %-13.3e|  %-13.3e|  %-13.3e|  %-13.3e|\n", p, z, w, t)
    end
end
dampreport(sys::LTISystem) = dampreport(STDOUT, sys)


@doc """`tzero(sys)`

Compute the invariant zeros of the system `sys`. If `sys` is a minimal
realization, these are also the transmission zeros.""" ->
function tzero(sys::TransferFunction)
    if issiso(sys)
        return roots(sys.matrix[1,1].num)
    else
        return tzero(ss(sys))
    end
end

# Implements the algorithm described in:
# Emami-Naeini, A. and P. Van Dooren, "Computation of Zeros of Linear
# Multivariable Systems," Automatica, 18 (1982), pp. 415–430.
# 
# Note that this returns either Vector{Complex64} or Vector{Float64}
function tzero(sys::StateSpace)
    # Balance the system
    A, B, C = balance_statespace(sys.A, sys.B, sys.C)
    D = sys.D

    # Compute a good tolerance
    meps = 10*eps()*norm([A B; C D])
    A, B, C, D = reduce_sys(A, B, C, D, meps)
    if isempty(A)   return Float64[]    end
    A, B, C, D = reduce_sys(A', C', B', D', meps)
    if isempty(A)   return Float64[]    end

    # Compress cols of [C D] to [0 Df]
    mat = [C D]
    # To ensure type-stability, we have to annote the type here, as qrfact
    # returns many different types.
    W = full(qrfact(mat')[:Q], thin=false)::Matrix{Float64}
    W = fliplr(W)
    mat = mat*W
    if fastrank(mat', meps) > 0
        nf = size(A, 1)
        m = size(D, 2)
        Af = ([A B] * W)[1:nf, 1:nf]
        Bf = ([eye(nf) zeros(nf, m)] * W)[1:nf, 1:nf]
        zs = eig(Af, Bf)[1]
    else
        zs = Float64[]
    end
    return zs
end

# Implements REDUCE in the Emami-Naeini & Van Dooren paper. Returns transformed
# A, B, C, D matrices. These are empty if there are no zeros.
function reduce_sys(A::Matrix{Float64}, B::Matrix{Float64}, C::Matrix{Float64},
        D::Matrix{Float64}, meps::Float64)
    Cbar, Dbar = C, D
    while true
        # Compress rows of D
        U = full(qrfact(D)[:Q], thin=false)::Matrix{Float64}
        D = U'*D
        C = U'*C
        sigma = fastrank(D, meps)
        Cbar = C[1:sigma, :]
        Dbar = D[1:sigma, :]
        Ctilde = C[(1 + sigma):end, :]
        if sigma == size(D, 1)
            break
        end

        # Compress columns of Ctilde
        V = full(qrfact(Ctilde')[:Q], thin=false)::Matrix{Float64}
        V = fliplr(V)
        Sj = Ctilde*V
        rho = fastrank(Sj', meps)
        nu = size(Sj, 2) - rho

        if rho == 0
            break
        elseif nu == 0
            # System has no zeros, return empty matrices
            A = B = Cbar = Dbar = Float64[]
            break
        end
        # Update System
        n, m = size(B)
        Vm = [V zeros(n, m); zeros(m, n) eye(m)]
        if sigma > 0
            M = [A B; Cbar Dbar]
            Vs = [V' zeros(n, sigma) ; zeros(sigma, n) eye(sigma)]
        else
            M = [A B]
            Vs = V'
        end
        sigma, rho, nu
        M = Vs * M * Vm
        A = M[1:nu, 1:nu]
        B = M[1:nu, (nu + rho + 1):end]
        C = M[(nu + 1):end, 1:nu]
        D = M[(nu + 1):end,  (nu + rho + 1):end]
    end
    return A, B, Cbar, Dbar
end

# Determine the number of non-zero rows, with meps as a tolerance. For an
# upper-triangular matrix, this is a good proxy for determining the row-rank.
function fastrank(A::Matrix{Float64}, meps::Float64)
    n, m = size(A, 1, 2)
    if n*m == 0     return 0    end
    norms = Array(Float64, n)
    for i = 1:n
        norms[i] = norm(A[i, :])
    end
    mrank = sum(norms .> meps)
    return mrank
end
