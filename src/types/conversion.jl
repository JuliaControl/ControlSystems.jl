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
Base.convert{T<:Number}(::Type{<:TransferFunction}, b::T) = tf([b])
Base.convert{T<:Number}(::Type{<:TransferFunction{<:SisoRational}}, b::T) = tf(b)
Base.convert{T<:Number}(::Type{<:TransferFunction{<:SisoZpk}}, b::T) = zpk(b)
Base.convert{T<:Number}(::Type{<:TransferFunction{<:SisoRational{T}}}, b::Number) = tf(T(b))
Base.convert{T<:Number}(::Type{<:TransferFunction{<:SisoZpk{T}}}, b::Number) = zpk(T(b))

function Base.convert(::Type{StateSpace{T1, MT}}, b::T2) where {T1, MT, T2<:Number}
    A = convert(MT, fill(zero(T1), (0,0)))
    B = convert(MT, fill(zero(T1), (0,1)))
    C = convert(MT, fill(zero(T1), (1,0)))
    D = convert(MT, fill(convert(T1, b), (1,1)))
    Ts = 0.0
    return StateSpace{T1,MT}(A,B,C,D,Ts)
end


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
    Gnew_matrix = Matrix{S}(size(G))
    for i in eachindex(G.matrix)
        Gnew_matrix[i] = convert(S, G.matrix[i])
    end
    return TransferFunction{S}(Gnew_matrix, G.Ts)
end

function convert(::Type{StateSpace{T,MT}}, sys::StateSpace) where {T, MT}
    return StateSpace{T,MT}(convert(MT, sys.A), convert(MT, sys.B), convert(MT, sys.C), convert(MT, sys.D), sys.Ts)
end


function Base.convert(::Type{<:StateSpace}, G::TransferFunction)
    if !isproper(G)
        error("System is improper, a state-space representation is impossible")
    end

    # TODO : These are added due to scoped for blocks, but is a hack. This
    # could be much cleaner.
    T0 = numeric_type(G)
    T = typeof(one(T0) / one(T0))

    Ac = Bc = Cc = Dc = A = B = C = D = Array{T}(0, 0)
    for i=1:ninputs(G)
        for j=1:noutputs(G)
            a, b, c, d = siso_tf_to_ss(G.matrix[j, i])
            if j > 1
                # vcat
                Ac = blkdiag(Ac, a)
                Bc = vcat(Bc, b)
                Cc = blkdiag(Cc, c)
                Dc = vcat(Dc, d)
            else
                Ac, Bc, Cc, Dc = a, b, c, d
            end
        end
        if i > 1
            # hcat
            A = blkdiag(A, Ac)
            B = blkdiag(B, Bc)
            C = hcat(C, Cc)
            D = hcat(D, Dc)
        else
            A, B, C, D = Ac, Bc, Cc, Dc
        end
    end
    # A, B, C = balance_statespace(A, B, C)[1:3] NOTE: Use balance?
    return ss(A, B, C, D, G.Ts)
end

siso_tf_to_ss(f::SisoTf) = siso_tf_to_ss(convert(SisoRational, f))

# Conversion to statespace on controllable canonical form
function siso_tf_to_ss(f::SisoRational)
    T0 = numeric_type(f)
    T = typeof(one(T0)/one(T0))

    num0, den0 = Base.num(f), Base.den(f)
    # Normalize the numerator and denominator
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
        A = diagm(fill(one(T), N-1), 1)
        A[end, :] .= -reverse(den)[1:end-1]

        B = fill(zero(T), N, 1)
        B[end] = one(T)

        C = Matrix{T}(1, N)
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
function balance_statespace(A::AbstractMatrix{P}, B::AbstractMatrix{P}, C::AbstractMatrix{P}, perm::Bool=false) where P <: BlasNumber
    nx = size(A, 1)
    nu = size(B, 2)
    ny = size(C, 1)

    # Compute the transformation matrix
    mag_A = abs.(A)
    mag_B = max.(abs.(B), 0)
    mag_C = max.(abs.(C), 0)
    T = balance_transform(mag_A, mag_B, mag_C, perm)

    # Perform the transformation
    A = T*A/T
    B = T*B
    C = C/T

    return A, B, C, T
end

# Fallback mehod for systems with exotic matrices (i.e. TrackedArrays)
balance_statespace(A::AbstractMatrix, B::AbstractMatrix, C::AbstractMatrix, args...) = A,B,C,I

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
    T = [A B; C zeros(size(C*B))]
    size(T,1) < size(T,2) && (T = [T; zeros(size(T,2)-size(T,1),size(T,2))])
    size(T,1) > size(T,2) && (T = [T zeros(size(T,1),size(T,1)-size(T,2))])
    S = diag(balance(T, false)[1])
    Sx = S[1:nx]
    Sio = S[nx+1]
    # Compute permutation of x (if requested)
    pvec = perm ? balance(A, true)[2] * [1:nx;] : [1:nx;]
    # Compute the transformation matrix
    T = zeros(promote_type(R, Float64), nx, nx)
    T[pvec, :] = Sio * diagm(1./Sx)
    return T
end

balance_transform(sys::StateSpace, perm::Bool=false) = balance_transform(sys.A,sys.B,sys.C,perm)


@doc """`sys = ss2tf(s::StateSpace)`, ` sys = ss2tf(A, B, C, Ts = 0; inputnames = "", outputnames = "")`

Convert a `StateSpace` realization to a `TransferFunction`""" ->
function convert(::Type{TransferFunction}, sys::StateSpace)
    #T = promote_type(primitivetype(A), Float64) # Int SS is not well represented by Int tf, (should actually be possible..)
    T0 = numeric_type(sys)
    T = typeof(sqrt(one(T0))) # NOTE: Reasonable?
    matrix = Matrix{SisoRational{T}}(size(sys))

    A, B, C, D = ssdata(sys)

    # The following follows from the matrix inversion lemma:
    # det(X + uᵀv) = det(X)(1 + vᵀX⁻¹u), or
    # det((sI-A)+BC) = det(sI-A)(1 + C(si-A)⁻¹B)
    charpolyA = charpoly(A)
    for i=1:ninputs(sys), j=1:noutputs(sys)
        num = charpoly(A-B*C) - charpolyA + D[j, i]*charpolyA
        matrix[j, i] = SisoRational{T}(num, charpolyA)
    end
    TransferFunction{SisoRational{T}}(matrix, get_Ts(sys))
end





# TODO: Could perhaps be made more accurate. See: An accurate and efficient
# algorithm for the computation of the # characteristic polynomial of a general square matrix.
function charpoly(A::AbstractMatrix{<:Complex})
    Λ = eigvals(A)
    T = promote_type(eltype(A), Complex128)

    return prod(roots2complex_poly_factors(Λ)) # Compute the polynomial factors directly?
end
function charpoly(A::AbstractMatrix{<:Real})
    Λ = eigvals(A)
    T = promote_type(eltype(A), Float64)

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
