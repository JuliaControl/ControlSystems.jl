function Base.convert(::Type{StateSpace}, t::TransferFunction)
    if !isproper(t)
        error("System is improper, a state-space representation is impossible")
    end
    ny, nu = size(t)
    mat = t.matrix
    # TODO : These are added due to scoped for blocks, but is a hack. This
    # could be much cleaner.
    Ac = Bc = Cc = Dc = A = B = C = D = Array(Float64, 0, 0)
    for i=1:nu
        for j=1:ny
            a, b, c, d = siso_tf_to_ss(mat[j, i])
            if j > 1
                # vcat
                Ac = blkdiag(Ac, a)
                Bc = vcat(Bc, b)
                Cc = blkdiag(Cc, c)
                Dc = vcat(Dc, d)
            else
                Ac, Bc, Cc, Dc = (a, b, c, d)
            end
        end
        if i > 1
            # hcat
            A = blkdiag(A, Ac)
            B = blkdiag(B, Bc)
            C = hcat(C, Cc)
            D = hcat(D, Dc)
        else
            A, B, C, D = (Ac, Bc, Cc, Dc)
        end
    end
    A, B, C = balance_statespace(A, B, C)[1:3]
    return ss(A, B, C, D, t.Ts, inputnames=t.inputnames, outputnames=t.outputnames)
end

Base.convert(::Type{StateSpace}, t::Real) = ss(t)

Base.promote_rule{T<:SisoTf}(::Type{StateSpace}, ::Type{TransferFunction{T}}) = StateSpace
Base.promote_rule{T<:Real}(::Type{StateSpace}, ::Type{T}) = StateSpace

siso_tf_to_ss(t::SisoTf) = siso_tf_to_ss(convert(SisoRational, t))

function siso_tf_to_ss(t::SisoRational)
    t = normalize_tf(t)
    tnum0 = numvec(t)
    tden = denvec(t)
    len = length(tden)
    #Pad tnum0 with zeros
    tnum = [zeros(len - length(tnum0)); tnum0]
    d = Array(Float64, 1, 1)
    d[1] = tnum[1]

    if len==1 || tnum == [0.0]
        a = zeros(0, 0)
        b = zeros(0, 1)
        c = zeros(1, 0)
    else
        tden = tden[2:end]
        a = [-tden' ; eye(len - 2, len - 1)]
        b = eye(len - 1, 1)
        c = tnum[2:len]' - d * tden[:]'
    end
    return float64mat(a), float64mat(b), float64mat(c), d
end

function normalize_tf(t::SisoRational)
    num = numpoly(t)
    den = denpoly(t)
    d = denvec(t)[1]
    return SisoTf(num/d, den/d)
end


"""
`A, B, C, T = balance_statespace{S}(A::Matrix{S}, B::Matrix{S}, C::Matrix{S}, perm::Bool=false)`

`sys, T = balance_statespace(sys::StateSpace, perm::Bool=false)`

Computes a balancing transformation `T` that attempts to scale the system so
that the row and column norms of [T*A/T T*B; C/T 0] are approximately equal.
If `perm=true`, the states in `A` are allowed to be reordered.

This is not the same as finding a balanced realization with equal and diagonal observability and reachability gramians, see `balreal`
"""
function balance_statespace{S}(A::Matrix{S}, B::Matrix{S},
        C::Matrix{S}, perm::Bool=false)
    nx = size(A, 1)
    nu = size(B,2)
    ny = size(C,1)

    # Compute the transformation matrix
    mag_A = abs(A)
    mag_B = maximum([abs(B)  zeros(S, nx, 1)], 2)
    mag_C = maximum([abs(C); zeros(S, 1, nx)], 1)
    T = balance_transform(mag_A, mag_B, mag_C, perm)

    # Perform the transformation
    A = T*A/T
    B = T*B
    C = C/T

    return A, B, C, T
end

function balance_statespace(sys::StateSpace, perm::Bool=false)
    A, B, C, T = balance_statespace(sys.A,sys.B,sys.C, perm)
    return ss(A,B,C,sys.D), T
end

"""
`T = balance_transform{R}(A::Matrix{R}, B::Matrix{R}, C::Matrix{R}, perm::Bool=false)`

`T = balance_transform(sys::StateSpace, perm::Bool=false) = balance_transform(A,B,C,perm)`

Computes a balancing transformation `T` that attempts to scale the system so
that the row and column norms of [T*A/T T*B; C/T 0] are approximately equal.
If `perm=true`, the states in `A` are allowed to be reordered.

This is not the same as finding a balanced realization with equal and diagonal observability and reachability gramians, see `balreal`
See also `balance_statespace`, `balance`
"""
function balance_transform{R}(A::Matrix{R}, B::Matrix{R}, C::Matrix{R}, perm::Bool=false)
    nx = size(A, 1)
    # Compute a scaling of the system matrix M
    S = diag(balance([A B; C zeros(size(C*B))], false)[1])
    Sx = S[1:nx]
    Sio = S[nx+1]
    # Compute permutation of x (if requested)
    pvec = perm ? balance(A, true)[2] * [1:nx;] : [1:nx;]
    # Compute the transformation matrix
    T = zeros(R, nx, nx)
    T[pvec, :] = Sio * diagm(1./Sx)
    return T
end

balance_transform(sys::StateSpace, perm::Bool=false) = balance_transform(sys.A,sys.B,sys.C,perm)


@doc """`sys = ss2tf(s::StateSpace)`, ` sys = ss2tf(A, B, C, Ts = 0; inputnames = "", outputnames = "")`

Convert a `StateSpace` realization to a `TransferFunction`""" ->
function ss2tf(s::StateSpace)
    return ss2tf(s.A, s.B, s.C, s.D, s.Ts, inputnames=s.inputnames, outputnames=s.outputnames)
end

function ss2tf(A, B, C, D, Ts = 0; inputnames = "", outputnames = "")
    nu,ny = size(B,2),size(C,1)
    ubernum = Matrix{Vector}(ny,nu)
    uberden = Matrix{Vector}(ny,nu)
    for i = 1:nu, j=1:ny
        ubernum[j,i],uberden[j,i] = sisoss2tf(A, B[:,i], C[j,:], D[j,i])
    end
    tf(ubernum,uberden, Ts, inputnames=inputnames, outputnames=outputnames)
end

function sisoss2tf(A, B, C, D)
    charpolA = charpoly(A)
    numP = charpoly(A-B*C) - charpolA + D*charpolA
    denP = charpolA
    return numP[1:length(numP)], denP[1:length(denP)]
end

tf(sys::StateSpace) = ss2tf(sys)
zpk(sys::StateSpace) = zpk(ss2tf(sys))

function charpoly(A)
    λ = eigvals(A);
    p = reduce(*,ControlSystems.Poly([1.]), ControlSystems.Poly[ControlSystems.Poly([1, -λᵢ]) for λᵢ in λ]);
    if maximum(imag(p[:])./(1+abs(real(p[:])))) < sqrt(eps())
        for i = 1:length(p)
            p[i] = real(p[i])
        end
    else
        error("Characteristic polynomial should be real")
    end
    p
end
