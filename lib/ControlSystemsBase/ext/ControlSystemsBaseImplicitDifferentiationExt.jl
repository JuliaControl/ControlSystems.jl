module ControlSystemsBaseImplicitDifferentiationExt
using ControlSystemsBase
import ControlSystemsBase: kalman, lqr, are, ContinuousType, DiscreteType
using ForwardDiff
using ForwardDiff: Dual
using ImplicitDifferentiation
using ComponentArrays
using LinearAlgebra


function forward_arec(pars)
    (; A,B,Q,R) = pars
    # Q = reshape(Q0, size(A))
    ControlSystemsBase.are(Continuous, A, B, Q, R), 0
end

function conditions_arec(pars, X, noneed)
    (; A,B,Q,R) = pars
    C = A'X
    C .+= C'
    C .+= Q
    XB = X*B
    mul!(C, XB, R\XB', -1, 1)
end

const implicit_arec = ImplicitFunction(forward_arec, conditions_arec)


"""
    are(::Continuous, A, B, Q, R::AbstractMatrix{<:Dual}; kwargs)

To make the ARE solver work with dual numbers, make sure that the `R` matrix has the dual number element type.
"""
function are(::ContinuousType, A::AbstractMatrix, B, Q, R::AbstractMatrix{<:Dual}; kwargs...)
    pars = ComponentVector(; A,B,Q,R)
    X0, _ = implicit_arec(pars)
    X = X0 isa AbstractMatrix ? X0 : reshape(X0, size(A))
    X
end

"""
    lqr(::Continuous, A, B, Q, R::AbstractMatrix{<:Dual})

To make the LQR solver work with dual numbers, make sure that the `R` matrix has the dual number element type.
"""
function lqr(::ContinuousType, A, B, Q, R::AbstractMatrix{<:Dual})
    X = are(Continuous, A, B, Q, R)
    R\(B'X)
end

"""
    kalman(::Continuous, A, C, Q, R::AbstractMatrix{<:Dual})

To make the Kalman solver work with dual numbers, make sure that the `R` matrix has the dual number element type.
"""
function kalman(::ContinuousType, A, C, Q, R::AbstractMatrix{<:Dual})
    X = are(Continuous, A', C', Q, R)
    (R\(C*X))'
end



## Discrete


function forward_ared(pars)
    (; A,B,Q,R) = pars
    # Q = reshape(Q0, size(A))
    ControlSystemsBase.are(Discrete, A, B, Q, R), 0
end

function conditions_ared(pars, X, noneed)
    # A'XA - X - (A'XB+S)(R+B'XB)^(-1)(B'XA+S') + Q = 0
    # A'X*A - X - (A'X*B)*((R+B'X*B)\(B'X*A)) + Q
    (; A,B,Q,R) = pars
    AX = A'X
    C = AX*A
    C .+= Q .- X
    AXB = AX*B
    C .-= AXB*(((B'X*B) .+= R)\AXB')
    C
end

const implicit_ared = ImplicitFunction(forward_ared, conditions_ared)


"""
    are(::Discrete, A, B, Q, R::AbstractMatrix{<:Dual}; kwargs)

To make the ARE solver work with dual numbers, make sure that the `R` matrix has the dual number element type.
"""
function are(::DiscreteType, A::AbstractMatrix, B, Q, R::AbstractMatrix{<:Dual}; kwargs...)
    pars = ComponentVector(; A,B,Q,R)
    X0, _ = implicit_ared(pars)
    X = X0 isa AbstractMatrix ? X0 : reshape(X0, size(A))
    X
end

"""
    lqr(::Discrete, A, B, Q, R::AbstractMatrix{<:Dual})

To make the LQR solver work with dual numbers, make sure that the `R` matrix has the dual number element type.
"""
function lqr(::DiscreteType, A, B, Q, R::AbstractMatrix{<:Dual})
    X = are(Discrete, A, B, Q, R)
    BX = B'X
    (R+BX*B)\(BX*A)
end

"""
    kalman(::Discrete, A, C, Q, R::AbstractMatrix{<:Dual})

To make the Kalman solver work with dual numbers, make sure that the `R` matrix has the dual number element type.
"""
function kalman(::DiscreteType, A, C, Q, R::AbstractMatrix{<:Dual})
    X = are(Discrete, A', C', Q, R)
    CX = C*X
    (R+CX*B)\(CX*A)
end


## Hinf norm
import ControlSystemsBase: hinfnorm
function forward_hinfnorm(pars; kwargs...)
    (; A,B,C,D) = pars
    sys = ss(A,B,C,D)
    hinfnorm(sys; kwargs...)
end

function conditions_hinfnorm(pars, γ, w; tol=1e-10)
    (; A,B,C,D) = pars
    sys = ss(A,B,C,D)
    [opnorm(freqresp(sys, w)) - γ]
end

const implicit_hinfnorm = ImplicitFunction(forward_hinfnorm, conditions_hinfnorm)

"""
    hinfnorm(sys::StateSpace{Continuous, <:Dual}; kwargs)

The H∞ norm can be differentiated through using ForwardDiff.jl, but at the time of writing, is limited to systems with *either* a signel input *or* a signle output. 

A reverse-differention rule is defined in RobustAndOptimalControl.jl, which means that hinfnorm is differentiable using, e.g., Zygote in reverse mode.
"""
function hinfnorm(sys::StateSpace{Continuous, <:Dual}; kwargs...)
    A,B,C,D = ssdata(sys)
    pars = ComponentVector(; A,B,C,D)
    γ, w = implicit_hinfnorm(pars)
    γ, w
end


# ## Schur, currently not working when the matrix A has complex eigenvalues.
# Sec 4.2 in "A PROCEDURE FOR DIFFERENTIATING PERFECT-FORESIGHT-MODEL REDUCED-FORM  OEFFICIENTS", Gary ANDERSON  has a formula for the derivative, but it looks rather expensive to compute, involving the factorization of a rather large Kronecker matrix. THis factorization only has to be done once, though, since it does not depend on the partials.
# function forward_schur(A)
#     F = schur(A)
#     ComponentVector(; F.Z, F.T), F
# end

# function conditions_schur(A, F, noneed)
#     (; Z, T) = F
#     [
#         vec(Z' * A * Z - T);
#         vec(Z' * Z - I + LowerTriangular(T) - Diagonal(T))
#     ]
# end

# linear_solver = (A, b) -> (Matrix(A) \ b, (solved=true,))
# const implicit_schur = ImplicitFunction(forward_schur, conditions_schur, linear_solver)

# # vectors = Z
# # Schur = T
# # A = F.vectors * F.Schur * F.vectors'
# # A = Z * T * Z'
# function LinearAlgebra.schur(A::AbstractMatrix{<:Dual})
#     ZT, F = implicit_schur(A)
#     n = length(A)
#     Z = reshape(ZT[1:n], size(A))
#     T = reshape(ZT[n+1:end], size(A))
#     Schur(T, Z, F.values)
# end


end # module

