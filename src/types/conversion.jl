# Base.convert(::Type{<:SisoTf}, b::Real) = Base.convert(SisoRational, b)
# Base.convert{T<:Real}(::Type{<:SisoZpk}, b::T) = SisoZpk(T[], T[], b)
# Base.convert{T<:Real}(::Type{<:SisoRational}, b::T) = SisoRational([b], [one(T)])
# Base.convert{T1}(::Type{SisoRational{Vector{T1}}}, t::SisoRational) =  SisoRational(Poly(T1.(t.num.a)),Poly(T1.(t.den.a)))
# Base.convert(::Type{<:StateSpace}, t::Real) = ss(t)
#

#
# function Base.convert{T<:AbstractMatrix{<:Number}}(::Type{StateSpace{T}}, s::StateSpace)
#     AT = promote_type(T, arraytype(s))
#     StateSpace{AT}(AT(s.A),AT(s.B),AT(s.C),AT(s.D), s.Ts, s.statenames, s.inputnames, s.outputnames)
# end

# TODO Fix these to use proper constructors
# NOTE: no real need to convert numbers to transfer functions, have addition methods..
# How to convert a number to either Continuous or Discrete transfer function
# In this case it would be motivated with a "Static" Type
# Base.convert(::Type{<:TransferFunction}, b::Number) = tf([b])
# Base.convert(::Type{<:TransferFunction{<:SisoRational}}, b::Number) = tf(b)
# Base.convert(::Type{<:TransferFunction{<:SisoZpk}}, b::Number) = zpk(b)

Base.convert(::Type{<:TransferFunction{<:SisoZpk{T1, TR1}}}, b::AbstractMatrix{T2}) where {T1, TR1, T2<:Number} = zpk(T1.(b))
Base.convert(::Type{<:TransferFunction{<:SisoRational{T1}}}, b::AbstractMatrix{T2}) where {T1, T2<:Number} = tf(T1.(b))
function convert(::Type{StateSpace{T,MT}}, D::AbstractMatrix{<:Number}) where {T, MT}
    (ny, nu) = size(D)
    A = MT(fill(zero(T), (0,0)))
    B = MT(fill(zero(T), (0,nu)))
    C = MT(fill(zero(T), (ny,0)))
    D = convert(MT, D)
    Ts = 0.0
    return StateSpace{T,MT}(A,B,C,D,Ts)
end

Base.convert(::Type{<:TransferFunction{<:SisoRational{T}}}, b::Number) where {T} = tf(T(b))
Base.convert(::Type{<:TransferFunction{<:SisoZpk{T, TR}}}, b::Number) where {T, TR} = zpk(T(b))
Base.convert(::Type{StateSpace{T, MT}}, b::Number) where {T, MT} = convert(StateSpace{T, MT}, fill(b, (1,1)))


#Base.convert(::Type{<:TransferFunction{<:SisoZpk}}, s::TransferFunction) = zpk(s)
#Base.convert(::Type{<:TransferFunction{<:SisoRational}}, s::TransferFunction) = tf(s)

#
# function Base.convert{T<:Real,S<:TransferFunction}(::Type{S}, b::VecOrMat{T})
#     r = Matrix{S}(size(b,2),1)
#     for j=1:size(b,2)
#         r[j] = vcat(map(k->convert(S,k),b[:,j])...)
#     end
#     hcat(r...)
# end
#

function convert(::Type{TransferFunction{S}}, G::TransferFunction) where S
    Gnew_matrix = convert.(S, G.matrix)
    return TransferFunction{eltype(Gnew_matrix)}(Gnew_matrix, G.Ts)
end

function convert(::Type{StateSpace{T,MT}}, sys::StateSpace) where {T, MT}
    return StateSpace{T,MT}(convert(MT, sys.A), convert(MT, sys.B), convert(MT, sys.C), convert(MT, sys.D), sys.Ts)
end

function Base.convert(::Type{StateSpace}, G::TransferFunction{<:SisoTf{T0}}) where {T0<:Number}
    T = Base.promote_op(/,T0,T0)
    convert(StateSpace{T,Matrix{T}}, G)
end


function Base.convert(::Type{StateSpace{T,MT}}, G::TransferFunction) where {T<:Number, MT<:AbstractArray{T}}
    if !isproper(G)
        error("System is improper, a state-space representation is impossible")
    end

    # TODO : These are added due to scoped for blocks, but is a hack. This
    # could be much cleaner.
    #T = Base.promote_op(/, T0, T0)

    Ac = Bc = Cc = Dc = A = B = C = D = Array{T}(undef, 0, 0)
    for i=1:ninputs(G)
        for j=1:noutputs(G)
            a, b, c, d = siso_tf_to_ss(T, G.matrix[j, i])
            if j > 1
                # vcat
                Ac = blockdiag(Ac, a)
                Bc = vcat(Bc, b)
                Cc = blockdiag(Cc, c)
                Dc = vcat(Dc, d)
            else
                Ac, Bc, Cc, Dc = a, b, c, d
            end
        end
        if i > 1
            # hcat
            A = blockdiag(A, Ac)
            B = blockdiag(B, Bc)
            C = hcat(C, Cc)
            D = hcat(D, Dc)
        else
            A, B, C, D = Ac, Bc, Cc, Dc
        end
    end
    # A, B, C = balance_statespace(A, B, C)[1:3] NOTE: Use balance?
    return StateSpace{T,MT}(A, B, C, D, G.Ts)
end

siso_tf_to_ss(T::Type, f::SisoTf) = siso_tf_to_ss(T, convert(SisoRational, f))

# Conversion to statespace on controllable canonical form
function siso_tf_to_ss(T::Type, f::SisoRational)

    num0, den0 = numvec(f), denvec(f)
    # Normalize the numerator and denominator,
    # To allow realization of transfer functions that are proper, but now strictly proper
    num = num0 / den0[1]
    den = den0 / den0[1]

    N = length(den) - 1 # The order of the rational function f

    # Get numerator coefficient of the same order as the denominator
    bN = length(num) == N+1 ? num[1] : 0

    if N==0 #|| num == zero(Poly{T})
        A = fill(zero(T), 0, 0)
        B = fill(zero(T), 0, 1)
        C = fill(zero(T), 1, 0)
    else
        A = diagm(1 => fill(one(T), N-1))
        A[end, :] .= -reverse(den)[1:end-1]

        B = fill(zero(T), N, 1)
        B[end] = one(T)

        C = fill(zero(T), 1, N)
        C[1:min(N, length(num))] = reverse(num)[1:min(N, length(num))]
        C[:] -= bN * reverse(den)[1:end-1] # Can index into polynomials at greater inddices than their length
    end
    D = fill(bN, 1, 1)

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
function balance_statespace(A::AbstractMatrix{P}, B::AbstractMatrix{P}, C::AbstractMatrix{P}, perm::Bool=false) where P <: BlasFloat
    nx = size(A, 1)
    nu = size(B, 2)
    ny = size(C, 1)

    # Compute the transformation matrix
    mag_A = abs.(A)
    mag_B = max.(abs.(B), zero(real(P)))
    mag_C = max.(abs.(C), zero(real(P)))
    T = balance_transform(mag_A, mag_B, mag_C, perm)

    # Perform the transformation
    A = T*A/T
    B = T*B
    C = C/T

    return A, B, C, T
end

#TODO Throw a warning here at least?
# Fallback mehod for systems with exotic matrices (i.e. TrackedArrays)
balance_statespace(A::AbstractMatrix, B::AbstractMatrix, C::AbstractMatrix, perm::Bool=false) = A,B,C,I

# First try to promote and see of we get BlasFloat, otherwise, fall beack on function above
function balance_statespace(A::Matrix{<:Number}, B::Matrix{<:Number}, C::Matrix{<:Number}, D::Matrix{<:Number}, perm::Bool=false)
    T = promote_type(eltype(A), eltype(B), eltype(C), eltype(D))
    A2, B2, C2, D2 = promote(A,B,C,D, fill(zero(T)/one(T),0,0)) # If Int, we get Float64
    balance_statespace(A2, B2, C2, perm)
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

No balancing will be done if `A, B C` are not BLAS compatible

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


convert(::Type{TransferFunction}, sys::StateSpace) = convert(TransferFunction{SisoRational}, sys)

function convert(::Type{TransferFunction{SisoRational{T}}}, sys::StateSpace) where {T<:Number}
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
    TransferFunction{SisoRational{T}}(matrix, get_Ts(sys))
end
function convert(::Type{TransferFunction{SisoRational}}, sys::StateSpace{T0}) where {T0<:Number}
    T = typeof(one(T0)/one(T0))
    convert(TransferFunction{SisoRational{T}}, sys)
end


function convert(::Type{TransferFunction{SisoZpk{T,TR}}}, sys::StateSpace) where {T<:Number, TR <: Number}
    matrix = Matrix{SisoZpk{T,TR}}(undef, size(sys))

    A, B, C, D = ssdata(sys)

    for i=1:noutputs(sys), j=1:ninputs(sys)
        z, p, k = siso_ss_to_zpk(sys, i, j)
        matrix[i, j] = SisoZpk{T,TR}(z, p, k)
    end
    TransferFunction{SisoZpk{T,TR}}(matrix, get_Ts(sys))
end
function convert(::Type{TransferFunction{SisoZpk}}, sys::StateSpace{T0}) where {T0<:Number}
    T = typeof(one(T0)/one(T0))
    convert(TransferFunction{SisoZpk{T,complex(T)}}, sys)
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
    Λ = eigvals(A)

    return prod(roots2poly_factors(Λ)) # Compute the polynomial factors directly?
end
function charpoly(A::AbstractMatrix{<:Real})
    Λ = eigvals(A)
    return prod(roots2real_poly_factors(Λ))
end


# function charpoly(A)
#     λ = eigvals(A);
#     T = promote_type(primitivetype(A), Float64)
#     I = one(T)
#     p = reduce(*,Poly([I]), Poly[Poly([I, -λᵢ]) for λᵢ in λ]);
#     if maximum(imag.(p[:])./(I+abs.(real.(p[:])))) < sqrt(eps(T))
#         for i = 1:length(p)
#             p[i] = real(p[i])
#         end
#     else
#         error("Characteristic polynomial should be real")
#     end
#     p
# end
