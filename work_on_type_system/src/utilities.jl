numeric_type(sys::StateSpace{T}) where T = T
numeric_type(sys::TransferFunction{<:SisoRational{T}}) where T = T
numeric_type(sys::TransferFunction{<:SisoZpk{T}}) where T = T


to_matrix(T, A::Vector) = Matrix{T}(reshape(A, length(A), 1))
to_matrix(T, A::AbstractMatrix) = Matrix{T}(A)
to_matrix(T, A::Number) = fill(T(A), 1, 1)


# NOTE: Tolerances for checking real-ness removed, shouldn't happen from LAPACK?
# TODO: This doesn't play too well with dual numbers..
# Allocate for maxiumum possible length of polynomial vector?
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
                error("Found pole without matching conjugate.")
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


"""
outs = index2range(ind1, ind2)
Helper function to convert indexes with scalars to ranges. Used to avoid dropping dimensions
"""
index2range(ind1, ind2) = (index2range(ind1), index2range(ind2))
index2range{T<:Number}(ind::T) = ind:ind
index2range{T<:AbstractArray}(ind::T) = ind
index2range(ind::Colon) = ind
