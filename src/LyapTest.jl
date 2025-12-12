module LyapTest

#=
The implementaiton is intended to be simple and fast. However in order to maintain
 readability, a few small optimizations have not been persued, for example:
* Symmetry is not exploited when solving diagonal entries of Lyapuniv equations
* Could use one more in-place multiplication in the top level functions
* Small savings could be possible when compting U*Q*U' where Q is symmetric, see [1]
* Could make the traiangular solve methods accepts various combinations of upper/lower
  Schur matrices but this seems to add quite a bit of code at quite a bit of code for quite a small gain.

# Some other questions
* What convention for signs and tranposes?
* Should error checking be done after each block solve or at the end of the algorithm?
=#

using LinearAlgebra
using StaticArrays



"""
    `_schurstructure(R::AbstractMatrix, ul=Union{Val{:U}, Val{:L}}) -> (b, d, nblocks)`

Return the block strucutre of an upper quasi-traingular Schur matrix `R`.
`ul` indicates if R is upper (`ul=Val(:U)`) or lower (`ul=Val(:L)`) triangular.


`d` contains the block sizees of each diagonal block (`1` or `2`)

`b` contains the indices of the blocks

`nblocks` is the number of blocks

"""
function _schurstructure(R::AbstractMatrix, ul=Val(:U)::Union{Val{:U}, Val{:L}})
    n = size(R,1)

    d = Vector{Int}(undef, n) # block sizes
    b = Vector{UnitRange{Int64}}(undef, n) # block indices

    j = 1 # column if ul=:U, row if ul=:L
    k = 0 # block number
    while j <= n
        k += 1
        if j == n
            d[k] = 1
        else
            if ul === Val(:U)
                d[k] = iszero(R[j+1, j]) ? 1 : 2
            else
                d[k] = iszero(R[j, j+1]) ? 1 : 2
            end
        end
        b[k] = j:j+d[k]-1
        j += d[k]
    end
    resize!(d, k)
    resize!(b, k)
    return d, b, k
end

# FIXME: better handling of uniform scaling?!
issquare(A::Number) = true
issquare(A::AbstractMatrix) = size(A,1) == size(A,2)
function _check_lyap_inputs(A, Q)
    if !issquare(A); error("The A matrix must be square"); end
    if !ishermitian(Q); error("The Q matrix must be Hermitian"); end
    if size(Q, 1) != size(A, 1); error("The A and Q matrices must have the same dimensions"); end
end

function _check_sylv_inputs(A, B, C)
    if !issquare(A); error("The A matrix must be square"); end
    if !issquare(B); error("The B matrix must be square"); end
    if size(C, 1) != size(A, 1); error("The A and C matrices have inconsistent dimensions"); end
    if size(C, 2) != size(B, 2); error("The B and C matrices have inconsistent dimensions"); end
end

# Should preferably be fixed in LinearAlgebra
LinearAlgebra.schur(A::AbstractMatrix{T}) where T = schur!(LinearAlgebra.copy_oftype(A, LinearAlgebra.eigtype(T)))


sylvcsoltype(A, B, C) = Base.promote_op((a,b,c) -> c / (a + b), eltype(A), eltype(B), eltype(C))
#sylvcsoltype(A, B, C) = Base.promote_op(sylvc, eltype(A), eltype(B), eltype(C))
sylvdsoltype(A, B, C) = Base.promote_op(sylvd, eltype(A), eltype(B), eltype(C))

"""
    _sylvc!(A, B, C) -> X

Find the solution `X` to the continuous-time Sylvester equation

`AX + XB = C`

for small matrices (1x1, 1x2, 2x1, 2x2), overwriting the input `C`.
"""
@inline function _sylvc!(A::AbstractMatrix, B::AbstractMatrix, C::AbstractMatrix)
    M, N = size(C)
    if M == 2 && N == 2
        _sylvc!(A, B, C, Val(2), Val(2))
    elseif M == 2 && N == 1
        _sylvc!(A, B, C, Val(2), Val(1))
    elseif M == 1 && N == 2
        _sylvc!(A, B, C, Val(1), Val(2))
    elseif M == 1 && N == 1
        _sylvc!(A, B, C, Val(1), Val(1))
    else
        error("Matrix dimensionsins should not be greater than 2")
    end
    return C
end
@inline function _sylvc!(A::AbstractMatrix, B::AbstractMatrix, C::AbstractMatrix, ::Val{M}, ::Val{N}) where {T <: Number, M, N}
    As = SMatrix{M,M}(A)
    Bs = SMatrix{N,N}(B)
    Cvs = SMatrix{M,N}(C)[:] # vectorization of C

    Xv = lu(kron(SMatrix{N,N}(I), As) + kron(transpose(Bs), SMatrix{M,M}(I))) \ Cvs # using the vectorization identity vec(AXB) = kron(B'*A)*vec(X) (with A = I or B = I)

    if any(!isfinite, Xv); error("Matrix equation has no solution, see ?sylvc or ?lyapc"); end

    C .= reshape(Xv, M, N)
end


"""
    _sylvd!(A, B, C) -> X

Find the solution `X` to the discrete-time Sylvester equation

`AXB - X = C`

for small matrices (1x1, 1x2, 2x1, 2x2), overwriting the input `C`.
"""
function _sylvd!(A::AbstractMatrix, B::AbstractMatrix, C::AbstractMatrix)
    M, N = size(C)
    if M == 2 && N == 2
        _sylvd!(A, B, C, Val(2), Val(2))
    elseif M == 2 && N == 1
        _sylvd!(A, B, C, Val(2), Val(1))
    elseif M == 1 && N == 2
        _sylvd!(A, B, C, Val(1), Val(2))
    elseif M == 1 && N == 1
        _sylvd!(A, B, C, Val(1), Val(1))
    else
        error("Matrix dimensionsins should not be greater than 2")
    end
    return C
end
function _sylvd!(A::AbstractMatrix, B::AbstractMatrix, C::AbstractMatrix, ::Val{M}, ::Val{N}) where {T <: Number, M, N}
    As = SMatrix{M,M}(A)
    Bs = SMatrix{N,N}(B)
    Cvs = SMatrix{M,N}(C)[:] # vectorization of C

    Xv = (kron(transpose(Bs), As) - SMatrix{M*N,M*N}(I)) \ Cvs # using the vectorization identity vec(AXB) = kron(B'*A)*vec(X)

    if any(!isfinite, Xv); error("Matrix equation has no solution, see ?sylvd or ?lyapd"); end

    C .= reshape(Xv, M, N)
end

"""
    sylvc(A, B, C) -> X

Find the solution `X` to the continuous-time Sylvester equation

`AX + XB = C`

A solution exists unless `A` and `B` have eigenvalues `λ` and `μ` such that λ + μ = 0.

[1] **Bartels, R. H., & Stewart, G. W.** (1972). "Solution of the matrix
    equation AX + XB = C" Communications of the ACM, 15(9), 820-826.
"""
function sylvc(A, B, C)
    _check_sylv_inputs(A, B, C)

    At2, UA = schur(A')
    B2, UB = schur(B)

    C2 = UA'*C*UB # This should give C2 the right type

    Y = _sylvc_schur!(Matrix(At2'), B2, C2, Val(:sylv))

    X = mul!(Y, UA, Y*UB')
end
@inline function sylvc(a::Number, b::Number, c::Number)
    x = c / (a + b)

    if !isfinite(x); error("Matrix equation has no solution, see ?sylvc or ?lyapc"); end

    x
end

"""
    sylvd(A, B, C) -> X

Find the solution `X` to the discrete-time Sylvester equation

`AXB - X = C`

A solution exists unless `A` and `B` have eigenvalues `λ` and `μ` such that λμ = 1.

[1] **Bartels, R. H., & Stewart, G. W.** (1972). "Solution of the matrix
    equation AX + XB = C" Communications of the ACM, 15(9), 820-826.
"""
function sylvd(A, B, C)
    _check_sylv_inputs(A, B, C)

    At2, UA = schur(A')
    B2, UB = schur(B)

    C2 = UA'*C*UB

    Y = _sylvd_schur!(Matrix(At2'), B2, C2, Val(:sylv))

    X = mul!(Y, UA, Y*UB')
end
@inline function sylvd(a::Number, b::Number, c::Number)
    x = c / (a * b - 1)

    if !isfinite(x); error("Matrix equation has no solution, see ?sylvd or ?lyapd"); end

    x
end

"""
    lyapc(A, Q) -> X

Computes the solution `X` of the continuous-time Lyapunov equation

`AX + XA' + Q = 0`

A solution exists unless `A` has an eigenvalue λ = ±1 or an eigenvalue pair λ₁λ₂ = 1.

[1] **Bartels, R. H., & Stewart, G. W.** (1972). "Solution of the matrix
    equation AX + XB = C" Communications of the ACM, 15(9), 820-826.
"""
function lyapc(A, Q)

     _check_lyap_inputs(A, Q)

    if !hasmethod(schur!, (typeof(A),))
        return lyapc(A, Q, Val(:naive))
    end

    At2, U = schur(A')

    Q2 = U'*Q*U

    Y = _sylvc_schur!(Matrix(At2'), At2, lmul!(-1, Q2), Val(:lyap))

    X = mul!(Y, U, Y*U')
end


function sylvc(A, B, C, ::Val{:naive})
    Xv = kron(I(size(B,1)), A) + kron(transpose(B), I(size(A,1))) \ C[:]
    return reshape(Xv, size(C))
end

# Mapping from Cartesian index into vector represenentation of Symmetric matrix
@inline sub2triidx(i,j) = (j >= i ? (j*(j-1))>>>1 + i : (i*(i-1))>>>1 + j)

function lyapc(A, Q, ::Val{:naive}) # Only works for real matrices A
    # Sets up and solves a system of equations for the upper triangular part of X
    # and solves that equation. This gives an n(n+1)/2 system instead of n^2
    # ONLY WORKS FOR REAL A!
    # Should be able to to base the discrete time version on this as well
    _check_lyap_inputs(A, Q)

    if !isreal(A); error("Only works for real A matrices"); end # Should call sylvc

    n = size(Q, 1)
    nR = (n*(n+1)) >>> 1 # Ssize of the system to be solved

    R = zeros(eltype(A), nR, nR)
    minusQv = Vector{eltype(Q)}(undef, nR)

    for j=1:n
        for i=1:j
            m = sub2triidx(i,j) # Set up equation for X[i,j]
            minusQv[m] = -Q[i,j]
            # (AX + XA')[i,j] = sum(A[i,k]*X[k,j]) + sum(X[i,k]*A'[k,j])
            # the r = trinum[i]+r gives the kth element in the upper traingle
            # which correpsonds to X[i,j]
            for k=1:n
                R[m, sub2triidx(k,j)] += A[i,k]
                R[m, sub2triidx(i,k)] += A[j,k]
            end
        end
    end

    Xv = R \ minusQv

    # Fill the upper traingle of X from Xv
    X = [Xv[sub2triidx(i,j)] for i=1:n, j=1:n]

    return X
end





"""
    X = lyapd(A, Q) -> X

Find the solution `X` to the discrete-time Lyapunov equation

`AXA' - X + Q = 0`

A solution exists unless `A` has an eigenvalue λ = ±1 or an eigenvalue pair λ₁λ₂ = 1 .


[1] **Barraud, A.** (1977) "A numerical algorithm to solve A'XA - X = Q"
    IEEE Transactions on Automatic Control

[2] **Bartels, R. H., & Stewart, G. W.** (1972). "Solution of the matrix
    equation AX + XB = C" Communications of the ACM, 15(9), 820-826.

"""
function lyapd(A, Q)

    _check_lyap_inputs(A, Q)

    At2, U = schur(A')

    Q2 = U'*Q*U

    Y = _sylvd_schur!(Matrix(At2'), At2, lmul!(-1, Q2), Val(:lyap))

    X = mul!(Y, U, Y*U')
end


"""
    sylvc_schur!(A::AbstractMatrix, B::AbstractMatrix, C::AbstractMatrix) -> X

Find the solution `X` to the continuous-time Sylvester equation

`AX + XB = C`

where `A` is assumed to have lower Schur form (quasi-triangular, 1x1 & 2x2 blocks on the diagonal)
`B` is assumed to have upper Schur form

See also `sylvc`
"""
function sylvc_schur!(A::AbstractMatrix, B::AbstractMatrix, C::AbstractMatrix)
    _check_sylv_inputs(A, B, C)
    T = sylvcsoltype(A, B, C)
    _sylvc_schur!(convert(Matrix, A), convert(Matrix, B), convert(Matrix{T}, C), Val(:sylv))
end
function lyapc_schur!(A::AbstractMatrix, Q::AbstractMatrix)
    _check_lyap_inputs(A, Q)
    T = sylvcsoltype(A, A, Q)
    _sylvc_schur!(convert(Matrix, A), convert(Matrix, A'), lmul!(-1, Matrix{T}(Q)), Val(:lyap))
end
function _sylvc_schur!(A::Matrix, B::Matrix, C::Matrix, alg::Union{Val{:sylv},Val{:lyap}},
    schurtype::Union{Val{:real},Val{:complex}} = isreal(A) || isreal(B) ? Val(:real) : Val(:complex)) where {T <: Number}
    # The user should preferably use sylvc_schur! and lyapc_schur!
    # I.e., this method does not check whether C is hermitian
    # The matrix C is successively replaced with the solution X
    # if alg === Val(:lyap), only the lower triangle of C is computed
    # after which an Hermitian view is applied

    # get block indices and nbr of blocks
    if schurtype === Val(:real)
        _, ba, nblocksa = _schurstructure(A, Val(:L)) # A is assumed upper triangualar
        _, bb, nblocksb = _schurstructure(B, Val(:U))
    else
        nblocksa = size(A, 1)
        nblocksb = size(B, 1)
    end

    @inbounds for j=1:nblocksb
        i0 = (alg === Val(:lyap) ? j : 1)

        # if j > 1; mul!(C[i0:nblocksa,j], C[i0:nblocksa, 1:j-1], B[1:j-1, j], 1, -1); end # Could move this out?
        # figure out the row indexing

        for i=i0:nblocksa
            if schurtype === Val(:complex)
                if i > 1; C[i,j] -= sum(A[i, k] * C[k, j] for k=1:i-1); end
                if j > 1; C[i,j] -= sum(C[i, k] * B[k, j] for k=1:j-1); end

                C[i,j] = sylvc(A[i, i], B[j, j], C[i, j]) # C[i,j] now contains  solution Y[i,j]

                if alg === Val(:lyap) && i > j
                    C[j,i] = conj(C[i,j])
                end
            else
                Aii = view(A, ba[i], ba[i])
                Bjj = view(B, bb[j], bb[j])
                Cij = view(C, ba[i], bb[j])

                if i > 1; @views mul!(Cij, A[ba[i], 1:ba[i-1][end]], C[1:ba[i-1][end], bb[j]], -1, 1); end
                if j > 1; @views mul!(Cij, C[ba[i], 1:bb[j-1][end]], B[1:bb[j-1][end], bb[j]], -1, 1); end

                _sylvc!(Aii, Bjj, Cij) # Cij now contains the solution Yij

                if alg === Val(:lyap) && i > j
                    for l=bb[j], k=ba[i]
                        C[l,k] = conj(C[k,l])
                    end
                end
            end
        end
    end
    return C
end


"""
    sylvd_schur!(A::AbstractMatrix, B::AbstractMatrix, C::AbstractMatrix) -> X

Solve the discrete-time Sylvester equation

`AXB - X = C`

where `A` is assumed to have lower Schur form (quasi-triangular, 1x1 & 2x2 blocks on the diagonal)
`B` is assumed to have upper Schur form

If the matrix `C` has the right type, it is overwritten with the solution `X`.

See also `sylvd`
"""
function sylvd_schur!(A::AbstractMatrix, B::AbstractMatrix, C::AbstractMatrix)
    _check_sylv_inputs(A, B, C)
    T = sylvdsoltype(A, B, C)
    _sylvd_schur!(convert(Matrix, A), convert(Matrix, B), convert(Matrix{T}, C), Val(:sylv))
end
function lyapd_schur!(A::AbstractMatrix, Q::AbstractMatrix)
    _check_lyap_inputs(A, Q)
    T = sylvdsoltype(A, A, Q)
    _sylvd_schur!(convert(Matrix, A), convert(Matrix, A'), lmul!(-1, Matrix{T}(Q)), Val(:lyap))
end
function _sylvd_schur!(A::Matrix, B::Matrix, C::Matrix, alg::Union{Val{:sylv},Val{:lyap}},
    schurtype::Union{Val{:real},Val{:complex}} = isreal(A) || isreal(B) ? Val(:real) : Val(:complex)) where {T <: Number}

    G = zeros(eltype(C), size(A,1), size(B, 1)) # Keep track of A*X for improved performance

    # get block dimensions, block indices, nbr of blocks
    if schurtype === Val(:real)
        _, ba, nblocksa = _schurstructure(A, Val(:L)) # A is assumed upper triangualar
        _, bb, nblocksb = _schurstructure(B, Val(:U))
    else
        nblocksa = size(A, 1)
        nblocksb = size(B, 1)
    end

    @inbounds for j=1:nblocksb
        i0 = (alg === Val(:lyap) ? j : 1)
        for i=i0:nblocksa
            if schurtype === Val(:complex)
                # Compute Gij up to the contribution from Aii*Yij which is added at the end of each iteration
                if i > 1; G[i,j] += sum(A[i,k] * C[k,j] for k=1:i-1); end

                C[i,j] -= sum(G[i,k] * B[k,j] for k=1:j)

                C[i,j] = sylvd(A[i,i], B[j,j], C[i,j]) # C[i,j] now contains  solution Y[i,j]

                if alg === Val(:lyap) && i > j
                    C[j,i] = conj(C[i,j])
                end

                G[i,j] += A[i, i] * C[i, j]
            else
                Aii = view(A, ba[i], ba[i])
                Bjj = view(B, bb[j], bb[j])
                Cij = view(C, ba[i], bb[j])

                Gij = view(G, ba[i], bb[j])

                if i > 1
                    @views mul!(Gij, A[ba[i], 1:ba[i-1][end]], C[1:ba[i-1][end], bb[j]], 1, 1)
                end

                @views mul!(Cij, G[ba[i], 1:bb[j][end]], B[1:bb[j][end], bb[j]], -1, 1)

                _sylvd!(Aii, Bjj, Cij) # Cij now contains the solution Yij

                if alg === Val(:lyap) && i > j
                    for l=bb[j], k=ba[i] # Avoids aliasing of copyto!
                        C[l,k] = conj(C[k,l])
                    end
                end

                mul!(Gij, Aii, Cij, 1, 1)
            end
        end
    end
    return C
end



end
