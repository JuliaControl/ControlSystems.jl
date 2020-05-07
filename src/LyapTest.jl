module LyapTest

using LinearAlgebra
using StaticArrays

#=
* Symmetry is not exploited when solving diagonal entries
* Perhaps it would be better with separate methods for the real and
  complex versions of the block solve algorithms
* What convention for signs and tranposes?
* Should error checking be done? In that case where?
=#

"""
    `b, d, nblocks = _schurstructure(R::AbstractMatrix)`

    Return the block strucutre of a quasi-traingular Schur matrix `R`.

    `d` contains the block sizees of each diagonal block (`1` or `2`)

    `b` contains the indices of the blocks

    `nblocks` is the number of blocks

"""
function _schurstructure(R::AbstractMatrix, ul=Val(:u)::Union{Val{:u}, Val{:l}})
    n = size(R,1)

    d = Vector{Int}(undef, n) # block sizes
    b = Vector{UnitRange{Int64}}(undef, n) # block indices

    j = 1 # column if ul=:u, row if ul=:l
    k = 0 # block number
    while j <= n
        k += 1
        if j == n
            d[k] = 1
        else
            if ul == Val(:u)
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


"""

    _sylvc!(A, B, C)

    Solve the continuous-time sylvester equation

    `AX + XB = C`

    for small matrices (1x1, 1x2, 2x1, 2x2)
"""
function _sylvc!(X::AbstractMatrix, A::SMatrix, B::SMatrix, C::SMatrix) where {T <: Number} # Should really mostly be used for Reals
    Cv = C[:] # vectorization of C
    m, n = size(C)
    Xv = (kron(I(n), A) + kron(transpose(B), I(m))) \ Cv # using the vectorization identity vec(AXB) = kron(B'*A)*vec(X)

    if any(!isfinite, Xv); error("Matrix equation has no solution, see ?sylvc or ?lyapc"); end

    X .= SMatrix{size(C,1),size(C,2)}(reshape(Xv, size(C,1), size(C,2)))
end
function _sylvc!(A::AbstractMatrix, B::AbstractMatrix, C::AbstractMatrix)
    if size(C) == (1,1)
        _sylvc!(C, SMatrix{1,1}(A), SMatrix{1,1}(B), SMatrix{1,1}(C))
    elseif size(C) == (1,2)
        _sylvc!(C, SMatrix{1,1}(A), SMatrix{2,2}(B), SMatrix{1,2}(C))
    elseif size(C) == (2,1)
        _sylvc!(C, SMatrix{2,2}(A), SMatrix{1,1}(B), SMatrix{2,1}(C))
    elseif size(C) == (2,2)
        _sylvc!(C, SMatrix{2,2}(A), SMatrix{2,2}(B), SMatrix{2,2}(C))
    else
        error("Matrix dimensionsins should not be greater than 2")
    end
    return C
end


"""

    _sylvd!(X, A, B, C)

    Solve the discrete-time sylvester equation

    `AXB - X = C`

    for small matrices (1x1, 1x2, 2x1, 2x2)
"""
function _sylvd!(A::AbstractMatrix, B::AbstractMatrix, C::AbstractMatrix)
    if size(C) == (1,1)
        _sylvd!(C, SMatrix{1,1}(A), SMatrix{1,1}(B), SMatrix{1,1}(C))
    elseif size(C) == (1,2)
        _sylvd!(C, SMatrix{1,1}(A), SMatrix{2,2}(B), SMatrix{1,2}(C))
    elseif size(C) == (2,1)
        _sylvd!(C, SMatrix{2,2}(A), SMatrix{1,1}(B), SMatrix{2,1}(C))
    elseif size(C) == (2,2)
        _sylvd!(C, SMatrix{2,2}(A), SMatrix{2,2}(B), SMatrix{2,2}(C))
    else
        error("Matrix dimensionsins should not be greater than 2")
    end
    return C
end
function _sylvd!(X::AbstractMatrix, A::SMatrix, B::SMatrix, C::SMatrix) where {T <: Number} # Should really mostly be used for Reals
    Cv = C[:] # vectorization of C
    Xv = (kron(transpose(B), A) - I) \ Cv # using the vectorization identity vec(AXB) = kron(B'*A)*vec(X)

    if any(!isfinite, Xv); error("Matrix equation has no solution, see ?sylvd or ?lyapd"); end

    X .= reshape(Xv, size(C,1), size(C,2))
end




"""
    Solve the continuous-time Sylvester equation

    `AX + XB = C`

    A solution `X` exists unless `A` and `B` have eigenvalues `λ` and `μ` such that λ + μ = 0.

    [1] **Bartels, R. H., & Stewart, G. W.** (1972). "Solution of the matrix
        equation AX + XB = C" Communications of the ACM, 15(9), 820-826.
"""
function sylvc(A, B, C)
    A2, UA = schur(Matrix(A'))
    B2, UB = schur(B)

    C2 = UA'*C*UB # C2 should have the right type

    Y = _sylvc_schur!(Matrix(A2'), B2, C2, Val(:sylv), isreal(A2) ? Val(:real) : Val(:complex))

    X = UA*Y*UB'
end
@inline function sylvc(a::Number, b::Number, c::Number)
    x = c / (a + b)

    if !isfinite(x); error("Matrix equation has no solution, see ?sylvc or ?lyapc"); end

    x
end

"""
    Solve the discrete-time Sylvester equation

    `AXB - X = C`

    A solution `X` exists unless `A` and `B` have eigenvalues `λ` and `μ` such that λμ = 1.

    [1] **Bartels, R. H., & Stewart, G. W.** (1972). "Solution of the matrix
        equation AX + XB = C" Communications of the ACM, 15(9), 820-826.
"""
function sylvd(A, B, C)
    A2, UA = schur(Matrix(A'))
    B2, UB = schur(B)

    C2 = UA'*C*UB

    Y = _sylvd_schur!(Matrix(A2'), B2, C2, Val(:sylv), isreal(A2) ? Val(:real) : Val(:complex))

    X = UA*Y*UB'
end
@inline function sylvd(a::Number, b::Number, c::Number)
    x = c / (a * b - 1)

    if !isfinite(x); error("Matrix equation has no solution, see ?sylvd or ?lyapd"); end

    x
end

"""
    Solve the continuous-time Lyapunov equation

    `AX + XA' + Q = 0`

    A solution exists if ... FIXME

    [1] **Bartels, R. H., & Stewart, G. W.** (1972). "Solution of the matrix
        equation AX + XB = C" Communications of the ACM, 15(9), 820-826.
"""
function lyapc(A, Q)

    if !ishermitian(Q); error("The Q matrix must be Hermitian"); end

    At2, U = schur(Matrix(A'))

    Q2 = U'*Q*U # Could use in-place, Small savings could be possible here, see [1]

    Y = _sylvc_schur!(Matrix(At2'), At2, -Q2, Val(:lyap), isreal(At2) ? Val(:real) : Val(:complex))

    X = U*Y*U' # Small savings could be possible here, see [1]
end
"""
    `X = lyapd(A, Q)`

    Find the solution `X` to the discrete-time Lyapunov equation

    `AXA' - X + Q = 0`

    A solution `X` exists if `A` has no eigenvalue λ = ±1 and no eigenvalue pair λ₁λ₂ = 1 .


    [1] **Barraud, A.** (1977) "A numerical algorithm to solve A'XA - X = Q"
        IEEE Transactions on Automatic Control

    [2] **Bartels, R. H., & Stewart, G. W.** (1972). "Solution of the matrix
        equation AX + XB = C" Communications of the ACM, 15(9), 820-826.

"""
function lyapd(A, Q)

    if !ishermitian(Q); error("The Q matrix must be Hermitian"); end

    At2, U = schur(Matrix(A'))

    Q2 = U'*Q*U # Some savings could be possible here, see [1]

    Y = _sylvd_schur!(Matrix(At2'), At2, -Q2, Val(:lyap), isreal(At2) ? Val(:real) : Val(:complex))

    X = U*Y*U' # Some savings could be possible here, see [1]
end


"""
    Solve the continuous-time Sylvetser equation

    `AX + XB = C`

    where
    `A` is assumed to have upper Schur form (lower quasi-triangular, 1x1 & 2x2 blocks on the diagonal)
    `B` is assumed to have lower Schur form

    If `alg == Val(:lyap)` then `C` should be Hermitian

    See also `sylvc`
"""
# It should also hold that `eltype(C) = eltype(C / (A + B))`
function _sylvc_schur!(A::AbstractMatrix, B::Matrix, C::AbstractMatrix, alg::Union{Val{:sylv},Val{:lyap}}, schurtype::Union{Val{:real},Val{:complex}}) where {T <: Number}
    # get block indices and nbr of blocks
    if schurtype == Val(:real)
        _, ba, nblocksa = _schurstructure(A, Val(:l)) # A is assumed upper triangualar
        _, bb, nblocksb = _schurstructure(B, Val(:u))
    else
        nblocksa = size(A, 1)
        nblocksb = size(B, 1)
    end

    for j=1:nblocksb
        i0 = (alg == Val(:lyap) ? j : 1)
        for i=i0:nblocksa
            if schurtype == Val(:complex)
                if i > 1; C[i,j] -= sum(A[i, k] * C[k, j] for k=1:i-1); end
                if j > 1; C[i,j] -= sum(C[i, k] * B[k, j] for k=1:j-1); end

                C[i,j] = sylvc(A[i, i], B[j, j], C[i, j]) # C[i,j] now contains  solution Y[i,j]

                if alg == Val(:lyap) && i > j
                    C[j,i] = conj(C[i,j])
                end
            else
                Cij = view(C, ba[i], bb[j])

                if i > 1; @views mul!(Cij, A[ba[i], 1:ba[i-1][end]], C[1:ba[i-1][end], bb[j]], -1, 1); end
                if j > 1; @views mul!(Cij, C[ba[i], 1:bb[j-1][end]], B[1:bb[j-1][end], bb[j]], -1, 1); end

                @views _sylvc!(A[ba[i], ba[i]], B[bb[j], bb[j]], Cij) # Cij now contains the solution Yij

                if alg == Val(:lyap) && i > j
                    view(C, ba[j], bb[i]) .= Cij'
                end
            end
        end
    end
    return C
end

# Should have eltype(C) = eltype(C / (A + B))
function _sylvd_schur!(A::AbstractMatrix, B::Matrix, C::AbstractMatrix, alg::Union{Val{:sylv},Val{:lyap}}, schurtype::Union{Val{:real},Val{:complex}}) where {T <: Number}
    # The matrix C is gradually replaced with the solution X

    G = zeros(eltype(C), size(A,1), size(B, 1)) # This matrix contains A*X for increased performance

    # get block dimensions, block indices, nbr of blocks
    if schurtype == Val(:real)
        _, ba, nblocksa = _schurstructure(A, Val(:l)) # A is assumed upper triangualar
        _, bb, nblocksb = _schurstructure(B, Val(:u))
    else
        nblocksa = size(A, 1)
        nblocksb = size(B, 1)
    end

    for j=1:nblocksb
        i0 = (alg == Val(:lyap) ? j : 1)
        for i=i0:nblocksa
            if schurtype == Val(:complex)
                # Compute Gij up to the contribution from Aii*Yij which is added at the end of each iteration
                if i > 1; G[i,j] += sum(A[i,k] * C[k,j] for k=1:i-1); end

                C[i,j] -= sum(G[i,k] * B[k,j] for k=1:j)

                C[i,j] = sylvd(A[i,i], B[j,j], C[i,j]) # C[i,j] now contains  solution Y[i,j]

                if alg == Val(:lyap) && i > j
                    C[j,i] = conj(C[i,j])
                end

                G[i,j] += A[i, i] * C[i, j]
            else
                Aii = A[ba[i], ba[i]]
                Bjj = B[bb[j], bb[j]]

                Cij = view(C, ba[i], bb[j])
                Gij = view(G, ba[i], bb[j])

                if i > 1
                    @views mul!(Gij, A[ba[i], 1:ba[i-1][end]], C[1:ba[i-1][end], bb[j]], 1, 1)
                end

                @views mul!(Cij, G[ba[i], 1:bb[j][end]], B[1:bb[j][end], bb[j]], -1, 1)

                _sylvd!(Aii, Bjj, Cij) # Cij now contains the solution Yij

                if alg == Val(:lyap) && i > j
                    view(C, ba[j], bb[i]) .= Cij'
                end

                mul!(Gij, Aii, Cij, 1, 1)
            end
        end
    end
    return C
end

end
