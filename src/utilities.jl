# Are only really needed for cases when we accept general LTISystem
# we should either use them consistently, with a good definition, or remove them
numeric_type(::Type{SisoRational{T}}) where T = T
numeric_type(::Type{<:SisoZpk{T}}) where T = T
numeric_type(sys::SisoTf) = numeric_type(typeof(sys))

numeric_type(::Type{TransferFunction{S}}) where S = numeric_type(S)
numeric_type(::Type{<:StateSpace{T}}) where T = T
numeric_type(sys::LTISystem) = numeric_type(typeof(sys))


to_matrix(T, A::Vector) = Matrix{T}(reshape(A, length(A), 1))
to_matrix(T, A::AbstractMatrix) = Matrix{T}(A)
to_matrix(T, A::Number) = fill(T(A), 1, 1)




# NOTE: Tolerances for checking real-ness removed, shouldn't happen from LAPACK?
# TODO: This doesn't play too well with dual numbers..
# Allocate for maxiumum possible length of polynomial vector?
#
# This function rely on that the every complex roots is followed by its exact conjugate,
# and that the first complex root in each pair has positive real part. This formaat is always
# returned by LAPACK routines for eigenvalues.
function roots2real_poly_factors(roots::Vector{cT}) where cT <: Number
    T = real(cT)
    poly_factors = Vector{Poly{T}}(0)

    for k=1:length(roots)
        r = roots[k]

        if isreal(r)
            push!(poly_factors,Poly{T}([-real(r),1]))
        else
            if imag(r) < 0 # This roots was handled in the previous iteration # TODO: Fix better error handling
                continue
            end

            if k == length(roots) || r != conj(roots[k+1])
                throw(AssertionError("Found pole without matching conjugate."))
            end

            push!(poly_factors,Poly{T}([real(r)^2+imag(r)^2, -2*real(r), 1]))
            # k += 1 # Skip one iteration in the loop
        end
    end

    return poly_factors
end

function roots2complex_poly_factors(roots::Vector{T}) where T <: Number
    return [Poly{T}([-r, 1]) for r in roots]
end

# Should probably try to get rif of this function...
poly2vec(p::Poly) = p.a[1:end]


function unwrap!(M::Array, dim=1)
    alldims(i) = [ n==dim ? i : (1:size(M,n)) for n in 1:ndims(M) ]
    for i = 2:size(M, dim)
        #This is a copy of slicedim from the JuliaLang but enables us to write to it
        #The code (with dim=1) is equivalent to
        # d = M[i,:,:,...,:] - M[i-1,:,...,:]
        # M[i,:,:,...,:] -= floor((d+π) / (2π)) * 2π
        d = M[alldims(i)...] - M[alldims(i-1)...]
        M[alldims(i)...] -= floor.((d+π) / 2π) * 2π
    end
    return M
end

#Collect will create a copy and collect the elements
unwrap(m::AbstractArray, args...) = unwrap!(collect(m), args...)
unwrap(x::Number) = x

"""
outs = index2range(ind1, ind2)
Helper function to convert indexes with scalars to ranges. Used to avoid dropping dimensions
"""
index2range(ind1, ind2) = (index2range(ind1), index2range(ind2))
index2range{T<:Number}(ind::T) = ind:ind
index2range{T<:AbstractArray}(ind::T) = ind
index2range(ind::Colon) = ind
