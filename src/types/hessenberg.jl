""" State Space where A is assumbed to be on Hessenberg form, only used from evalfr """
struct HessenbergSS{MT <: AbstractMatrix{<:Number}}
    A::MT
    B::MT
    C::MT
    D::MT
end

function evalfr(sys::HessenbergSS, s::Number)
    H = complex(sys.A)
    n = size(H,1)
    # Remove Iiw from H
    for i = 1:n
        H[i,i] -= im*w
    end
    rhs = complex(sys.B)
    backsolve_hess!(H, rhs)
    scale!(rhs, -1)
    return sys.D + sys.C*rhs
end


""" backsolve_hess!(H, rhs)
    Solve Hx=rhs for hessenberg H
"""
function backsolve_hess!(H::AbstractMatrix{T}, rhs::AbstractArray{T}) where {T}
    # Make H UpperTriangular and apply givens to rhs
    n = size(H,1)
    for i = 1:n-1
        Gi, r = Base.LinAlg.givens(H, i, i+1, i)
        A_mul_B!(Gi, H)
        A_mul_B!(Gi, rhs)
    end
    # H is now upper triangular
    A_ldiv_B!(UpperTriangular(H), rhs)
    return rhs
end
#
# ### Used with factorization of A if B and C are not changed
# """ solvehess(F, B, w)
#     Solves
#  (Iw-A)⁻¹B = Q(Iw-H)⁻¹Q⁻¹B
# where QHQ⁻¹ is a Hessenberg factoriztion
# """
# function solvehess(F::Base.LinAlg.Hessenberg{Tfact}, B, w) where {Tfact}
#     # For givens rotations Gₙ...G₂G₁ that turn Iw-H=V upper triangular
#     # Gₙ...G₂G₁V=U, we get
#     # QV⁻¹Q⁻¹B = Q(G₁ᵀG₂ᵀ...GₙᵀU)⁻¹Q⁻¹B = QU⁻¹Gₙ...G₂G₁Q⁻¹B
#     m,n = size(F.factors)
#     Q = F[:Q]::Base.LinAlg.HessenbergQ{Tfact,Array{Tfact,2}} #This call is type-unstable
#
#     T = typeof(im/one(Tfact)) # im is Complex{Bool}, T=Complex{Tfact} if Tfact::Real
#     H = Array{T, 2}(m,n)
#     copy!(H, F.factors)
#     triu!(H, -1) # H is now hessenberg and complex
#
#     # Remove Iiw from H
#     for i = 1:n
#         H[i,i] -= im*w
#     end
#
#     rhs = complex(Ac_mul_B(Q, B))
#
#     backsolve_hess!(H, rhs)
#
#     # Do Q*rhs
#     if Tfact == T
#         # Both Q and rhs should be complex
#         A_mul_B!(Q, rhs)
#     else
#         # Q is probably real (if A real)
#         # so do real multiplications to make it efficient
#         tmp1 = real(rhs)
#         tmp2 = imag(rhs)
#         A_mul_B!(Q, tmp1)
#         A_mul_B!(Q, tmp2)
#         rhs .= tmp1 .+ im.*tmp2
#     end
#
#     # We just solved for (A-Iiw), not (Iiw-A)
#     scale!(rhs, -1)
#     return rhs
# end
