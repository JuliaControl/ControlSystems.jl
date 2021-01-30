# Base.convert(::Type{<:SisoTf}, b::Real) = Base.convert(SisoRational, b)
# Base.convert{T<:Real}(::Type{<:SisoZpk}, b::T) = SisoZpk(T[], T[], b)
# Base.convert{T<:Real}(::Type{<:SisoRational}, b::T) = SisoRational([b], [one(T)])
# Base.convert{T1}(::Type{SisoRational{Vector{T1}}}, t::SisoRational) =  SisoRational(Polynomial(T1.(t.num.coeffs$1),Polynomial(T1.(t.den.coeffs$1))
# Base.convert(::Type{<:StateSpace}, t::Real) = ss(t)
#

#
# function Base.convert{T<:AbstractMatrix{<:Number}}(::Type{StateSpace{T}}, s::StateSpace)
#     AT = promote_type(T, arraytype(s))
#     StateSpace{AT}(AT(s.A),AT(s.B),AT(s.C),AT(s.D), s.timeevol, s.statenames, s.inputnames, s.outputnames)
# end

# TODO Fix these to use proper constructors
# NOTE: no real need to convert numbers to transfer functions, have addition methods..
# How to convert a number to either Continuous or Discrete transfer function
# In this case it would be motivated with a "Static" Type
# Base.convert(::Type{<:TransferFunction}, b::Number) = tf([b])
# Base.convert(::Type{<:TransferFunction{<:SisoRational}}, b::Number) = tf(b)
# Base.convert(::Type{<:TransferFunction{<:SisoZpk}}, b::Number) = zpk(b)
#
Base.convert(::Type{TransferFunction{TE,SisoZpk{T1, TR1}}}, D::AbstractMatrix{T2}) where {TE, T1, TR1, T2<:Number} =
    zpk(convert(Matrix{T1}, D), undef_sampletime(TE))
Base.convert(::Type{TransferFunction{TE,SisoRational{T1}}}, D::AbstractMatrix{T2}) where {TE, T1, T2<:Number} =
    tf(Matrix{T1}(D), undef_sampletime(TE))

function convert(::Type{StateSpace{TE,T}}, D::AbstractMatrix{<:Number}) where {TE,T}
    (ny, nu) = size(D)
    A = fill(zero(T), (0,0))
    B = fill(zero(T), (0,nu))
    C = fill(zero(T), (ny,0))
    D = convert(Matrix{T}, D)
    return StateSpace{TE,T}(A,B,C,D,undef_sampletime(TE))
end

# TODO We still dont have TransferFunction{..}(number/array) constructors
Base.convert(::Type{TransferFunction{TE,SisoRational{T}}}, d::Number) where {TE, T} =
    tf(T(d), undef_sampletime(TE))

Base.convert(::Type{TransferFunction{TE,SisoZpk{T,TR}}}, d::Number) where {TE, T, TR} =
    zpk(T(d), undef_sampletime(TE))

Base.convert(::Type{StateSpace{TE,T}}, d::Number) where {TE, T} =
    convert(StateSpace{TE,T}, fill(d, (1,1)))
#
# Base.convert(::Type{TransferFunction{Continuous,<:SisoRational{T}}}, b::Number) where {T} = tf(T(b), Continuous())
# Base.convert(::Type{TransferFunction{Continuous,<:SisoZpk{T, TR}}}, b::Number) where {T, TR} = zpk(T(b), Continuous())
# Base.convert(::Type{StateSpace{Continuous,T, MT}}, b::Number) where {T, MT} = convert(StateSpace{Continuous,T, MT}, fill(b, (1,1)))

# function Base.convert{T<:Real,S<:TransferFunction}(::Type{S}, b::VecOrMat{T})
#     r = Matrix{S}(size(b,2),1)
#     for j=1:size(b,2)
#         r[j] = vcat(map(k->convert(S,k),b[:,j])...)
#     end
#     hcat(r...)
# end
#

function convert(::Type{TransferFunction{TE,S}}, G::TransferFunction) where {TE,S}
    Gnew_matrix = convert.(S, G.matrix)
    return TransferFunction{TE,eltype(Gnew_matrix)}(Gnew_matrix, TE(G.timeevol))
end

function convert(::Type{StateSpace{TE,T}}, sys::AbstractStateSpace) where {TE, T}
    if sys isa StateSpace{TE, T}
        return sys
    else
        return StateSpace{TE, T}(convert(Matrix{T}, sys.A), convert(Matrix{T}, sys.B), convert(Matrix{T}, sys.C), convert(Matrix{T}, sys.D), TE(sys.timeevol))
    end
end

# TODO Maybe add convert on matrices
Base.convert(::Type{HeteroStateSpace{TE1,AT,BT,CT,DT}}, s::StateSpace{TE2,T}) where {TE1,TE2,T,AT,BT,CT,DT} =
    HeteroStateSpace{TE1,AT,BT,CT,DT}(s.A,s.B,s.C,s.D,TE1(s.timeevol))

Base.convert(::Type{HeteroStateSpace}, s::StateSpace) = HeteroStateSpace(s)

Base.convert(::Type{StateSpace}, s::HeteroStateSpace) = StateSpace(s.A, s.B, s.C, s.D, s.Ts)
Base.convert(::Type{StateSpace}, s::HeteroStateSpace{Continuous}) = StateSpace(s.A, s.B, s.C, s.D)

function Base.convert(::Type{StateSpace}, G::TransferFunction{TE,<:SisoTf{T0}}) where {TE,T0<:Number}

    T = Base.promote_op(/,T0,T0)
    convert(StateSpace{TE,T}, G)
end

function Base.convert(::Type{StateSpace{TE,T}}, G::TransferFunction) where {TE,T<:Number}
    if !isproper(G)
        error("System is improper, a state-space representation is impossible")
    end

    ny, nu = size(G)

    # A, B, C, D matrices for each element of the transfer function matrix
    abcd_vec = [siso_tf_to_ss(T, g) for g in G.matrix[:]]

    # Number of states for each transfer function element realization
    nvec = [size(abcd[1], 1) for abcd in abcd_vec]
    ntot = sum(nvec)

    A = zeros(T, (ntot, ntot))
    B = zeros(T, (ntot, nu))
    C = zeros(T, (ny, ntot))
    D = zeros(T, (ny, nu))

    inds = -1:0
    for j=1:nu
        for i=1:ny
            k = (j-1)*ny + i

            # states correpsonding to the transfer function element (i,j)
            inds = (inds.stop+1):(inds.stop+nvec[k])

            A[inds,inds], B[inds,j:j], C[i:i,inds], D[i:i,j:j] = abcd_vec[k]
        end
    end
    # A, B, C = balance_statespace(A, B, C)[1:3] NOTE: Use balance?
    return StateSpace{TE,T}(A, B, C, D, TE(G.timeevol))
end

siso_tf_to_ss(T::Type, f::SisoTf) = siso_tf_to_ss(T, convert(SisoRational, f))

# Conversion to statespace on controllable canonical form
function siso_tf_to_ss(T::Type, f::SisoRational)

    num0, den0 = numvec(f), denvec(f)
    # Normalize the numerator and denominator to allow realization of transfer functions
    # that are proper, but not strictly proper
    num = num0 / den0[1]
    den = den0 / den0[1]

    N = length(den) - 1 # The order of the rational function f

    # Get numerator coefficient of the same order as the denominator
    bN = length(num) == N+1 ? num[1] : 0

    if N == 0 #|| num == zero(Polynomial{T})
        A = zeros(T, (0, 0))
        B = zeros(T, (0, 1))
        C = zeros(T, (1, 0))
    else
        A = diagm(1 => ones(T, N-1))
        A[end, :] .= -reverse(den)[1:end-1]

        B = zeros(T, (N, 1))
        B[end] = one(T)

        C = zeros(T, (1, N))
        C[1:min(N, length(num))] = reverse(num)[1:min(N, length(num))]
        C[:] -= bN * reverse(den)[1:end-1] # Can index into polynomials at greater inddices than their length
    end
    D = fill(bN, (1, 1))

    return A, B, C, D
end

"""
`A, B, C, T = balance_statespace{S}(A::Matrix{S}, B::Matrix{S}, C::Matrix{S}, perm::Bool=false)`

`sys, T = balance_statespace(sys::StateSpace, perm::Bool=false)`

Computes a balancing transformation `T` that attempts to scale the system so
that the row and column norms of [T*A/T T*B; C/T 0] are approximately equal.
If `perm=true`, the states in `A` are allowed to be reordered.

This is not the same as finding a balanced realization with equal and diagonal observability and reachability gramians, see `balreal`
"""
function balance_statespace(A::AbstractMatrix, B::AbstractMatrix, C::AbstractMatrix, perm::Bool=false)
    try
        return _balance_statespace(A,B,C, perm)
    catch
        @warn "Unable to balance state-space, returning original system"
        return A,B,C,I
    end
 end

# # First try to promote and hopefully get some types we can work with
# function balance_statespace(A::AbstractMatrix, B::AbstractMatrix, C::AbstractMatrix, perm::Bool=false)
#     T = promote_type(eltype(A), eltype(B), eltype(C))
#     A2, B2, C2, D2 = promote(A,B,C, fill(zero(T)/one(T),0,0)) # If Int, we get Float64
#     balance_statespace(A2, B2, C2, perm)
# end

function balance_statespace(sys::StateSpace, perm::Bool=false)
    A, B, C, T = balance_statespace(sys.A,sys.B,sys.C, perm)
    return ss(A,B,C,sys.D,sys.timeevol), T
end

# Method that might fail for some exotic types, such as TrackedArrays
function _balance_statespace(A::AbstractMatrix, B::AbstractMatrix, C::AbstractMatrix, perm::Bool=false)
    nx = size(A, 1)
    nu = size(B, 2)
    ny = size(C, 1)

    # Compute the transformation matrix
    mag_A = abs.(A)
    mag_B = max.(abs.(B), false) # false is 0 of lowest type
    mag_C = max.(abs.(C), false)
    T = balance_transform(mag_A, mag_B, mag_C, perm)

    # Perform the transformation
    A = T*A/T
    B = T*B
    C = C/T

    return A, B, C, T
end

"""
`T = balance_transform{R}(A::AbstractArray, B::AbstractArray, C::AbstractArray, perm::Bool=false)`

`T = balance_transform(sys::StateSpace, perm::Bool=false) = balance_transform(A,B,C,perm)`

Computes a balancing transformation `T` that attempts to scale the system so
that the row and column norms of [T*A/T T*B; C/T 0] are approximately equal.
If `perm=true`, the states in `A` are allowed to be reordered.

This is not the same as finding a balanced realization with equal and diagonal observability and reachability gramians, see `balreal`
See also `balance_statespace`, `balance`
"""
function balance_transform(A::AbstractArray, B::AbstractArray, C::AbstractArray, perm::Bool=false)
    nx = size(A, 1)
    # Compute a scaling of the system matrix M
    R = promote_type(eltype(A), eltype(B), eltype(C), Float32) # Make sure we get at least BlasFloat
    T = R[A B; C zeros(R, size(C*B))]

    size(T,1) < size(T,2) && (T = [T; zeros(R, size(T,2)-size(T,1),size(T,2))])
    size(T,1) > size(T,2) && (T = [T zeros(R, size(T,1),size(T,1)-size(T,2))])
    S = diag(balance(T, false)[1])
    Sx = S[1:nx]
    Sio = S[nx+1]
    # Compute permutation of x (if requested)
    pvec = perm ? balance(A, true)[2] * [1:nx;] : [1:nx;]
    # Compute the transformation matrix
    T = zeros(R, nx, nx)
    T[pvec, :] = Sio * diagm(0 => R(1)./Sx)
    return T
end

balance_transform(sys::StateSpace, perm::Bool=false) = balance_transform(sys.A,sys.B,sys.C,perm)


convert(::Type{TransferFunction}, sys::AbstractStateSpace{TE}) where TE = convert(TransferFunction{TE,SisoRational}, sys)

function convert(::Type{TransferFunction{TE,SisoRational{T}}}, sys::AbstractStateSpace) where {TE,T<:Number}
    matrix = Matrix{SisoRational{T}}(undef, size(sys))

    A, B, C, D = ssdata(sys)

    # The following follows from the matrix inversion lemma:
    # det(X + uᵀv) = det(X)(1 + vᵀX⁻¹u), or
    # det((sI-A)+BC) = det(sI-A)(1 + C(si-A)⁻¹B)
    # C(si-A)⁻¹B) + D = 1/det(sI-A) * (det((sI-A)+BC) - I + D*det(sI-A))
    charpolyA = charpoly(A)
    for i=1:ninputs(sys), j=1:noutputs(sys)
        num = charpoly(A-B[:,i:i]*C[j:j,:]) - charpolyA + D[j, i]*charpolyA
        matrix[j, i] = SisoRational{T}(num, charpolyA)
    end
    TransferFunction{TE,SisoRational{T}}(matrix, TE(sys.timeevol))
end
function convert(::Type{TransferFunction{TE1,SisoRational}}, sys::StateSpace{TE2,T0}) where {TE1,TE2,T0<:Number}
    T = typeof(one(T0)/one(T0))
    convert(TransferFunction{TE1,SisoRational{T}}, sys)
end


function convert(::Type{TransferFunction{TE,SisoZpk{T,TR}}}, sys::StateSpace) where {TE,T<:Number, TR <: Number}
    matrix = Matrix{SisoZpk{T,TR}}(undef, size(sys))

    for i=1:noutputs(sys), j=1:ninputs(sys)
        z, p, k = siso_ss_to_zpk(sys, i, j)
        matrix[i, j] = SisoZpk{T,TR}(z, p, k)
    end
    TransferFunction{TE,SisoZpk{T,TR}}(matrix, TE(sys.timeevol))
end
function convert(::Type{TransferFunction{TE1,SisoZpk}}, sys::StateSpace{TE2,T0}) where {TE1,TE2,T0<:Number}
    T = typeof(one(T0)/one(T0))
    convert(TransferFunction{TE1,SisoZpk{T,complex(T)}}, sys)
end

"""
Convert get zpk representation of sys from input j to output i
"""
function siso_ss_to_zpk(sys, i, j)
    A, B, C = struct_ctrb_obsv(sys.A, sys.B[:, j:j], sys.C[i:i, :])
    D = sys.D[i:i, j:j]
    z = tzero(A, B, C, D)
    nx = size(A, 1)
    nz = length(z)
    k = nz == nx ? D[1] : (C*(A^(nx - nz - 1))*B)[1]
    return z, eigvals(A), k
end


# TODO: Could perhaps be made more accurate. See: An accurate and efficient
# algorithm for the computation of the # characteristic polynomial of a general square matrix.
function charpoly(A::AbstractMatrix{<:Number})
    Λ = eigvalsnosort(A)

    return prod(roots2poly_factors(Λ)) # Compute the polynomial factors directly?
end
function charpoly(A::AbstractMatrix{<:Real})
    Λ = eigvalsnosort(A)
    return prod(roots2real_poly_factors(Λ))
end


# function charpoly(A)
#     λ = eigvals(A);
#     T = promote_type(primitivetype(A), Float64)
#     I = one(T)
#     p = reduce(*,Polynomial([I]), Polynomial[Polynomial([I, -λᵢ]) for λᵢ in λ]);
#     if maximum(imag.(p[:])./(I+abs.(real.(p[:])))) < sqrt(eps(T))
#         for i = 1:length(p)
#             p[i] = real(p[i])
#         end
#     else
#         error("Characteristic polynomial should be real")
#     end
#     p
# end
