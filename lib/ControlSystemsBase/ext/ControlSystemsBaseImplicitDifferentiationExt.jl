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
    ControlSystemsBase.are(Continuous, A, B, Q, R), 0
end

function conditions_arec(pars, X, noneed)
    (; A,B,Q,R) = pars
    C = A'X
    C .+= C'
    C .+= Q
    XB = X*B
    mul!(C, XB, R\XB', -1, 1)
    # C .+ X .- X' # Does not seem to be needed
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
    C .+= Q .- X'
    AXB = AX*B
    C .-= AXB*(((B'X*B) .+= R)\AXB')
    # C .+= X .- X'
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
    (R+CX*C')\(CX*A')
end


## Lyap
function forward_lyapc(pars)
    (; A,Q) = pars
    # Q = reshape(Q0, size(A))
    ControlSystemsBase.lyapc(A, Q), 0
end

function conditions_lyapc(pars, X, noneed)
    (; A,Q) = pars
    AX = A*X
    O = AX .+ AX' .+ Q
    vec(O) + vec(X - X')
end

# linear_solver = (A, b) -> (Matrix(A) \ b, (solved=true,))
const implicit_lyapc = ImplicitFunction(forward_lyapc, conditions_lyapc)

"""
    ControlSystemsBase.lyap(nothing::ContinuousType, A::AbstractMatrix, Q::AbstractMatrix{<:Dual}; kwargs)

To make the Lyapunov solver work with dual numbers, make sure that the `Q` matrix has the dual number element type.

The returned gradient may not be symmetric, but the trick `(X + X') ./ 2` results in the correct symmetric gradient.

# Example:
```julia
using ControlSystemsBase, ImplicitDifferentiation, ForwardDiff, FiniteDifferences, ComponentArrays, Test
fdgrad(f, x) = FiniteDifferences.grad(central_fdm(3, 1), f, x) |> first

P = ssrand(1, 1, 2)
function difffun(q)
    Q = reshape(q, 2, 2)
    sum(ControlSystemsBase.lyap(P, Q))
end

q = [2.0 1; 1 2] |> vec
J1 = ForwardDiff.gradient(difffun, q) # Non-symmetric gradient

J1 = reshape(J1, 2,2)
J1 = vec((J1 + J1') ./ 2) # Symmetric gradient
J2 = fdgrad(difffun, q)
@test J1 ≈ J2 rtol = 1e-6
```
"""
function ControlSystemsBase.lyap(::ContinuousType, A::AbstractMatrix, Q::AbstractMatrix{<:Dual}; kwargs...)
    pars = ComponentVector(; A,Q)
    X0, _ = implicit_lyapc(pars)
    X0 isa AbstractMatrix ? X0 : reshape(X0, size(A))
end


# plyap
function forward_plyapc(pars)
    (; A,Q) = pars
    ControlSystemsBase.plyapc(A, Q), 0
end

function conditions_plyapc(pars, Xc, noneed)
    (; A,Q) = pars
    Q = Q*Q'
    X = Xc*Xc'
    AX = A*X
    O = AX .+ AX' .+ Q
    vec(O) + vec(Xc - UpperTriangular(Xc))
end

# linear_solver = (A, b) -> (Matrix(A) \ b, (solved=true,))
const implicit_plyapc = ImplicitFunction(forward_plyapc, conditions_plyapc)

function ControlSystemsBase.plyap(::ContinuousType, A::AbstractMatrix, Q::AbstractMatrix{<:Dual}; kwargs...)
    pars = ComponentVector(; A,Q)
    X0, _ = implicit_plyapc(pars)
    X0 isa AbstractMatrix ? X0 : reshape(X0, size(A))
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

The H∞ norm can be differentiated through using ForwardDiff.jl, but at the time of writing, is limited to systems with *either* a single input *or* a single output. 

A reverse-differentiation rule is defined in RobustAndOptimalControl.jl, which means that hinfnorm is differentiable using, e.g., Zygote in reverse mode.
"""
function hinfnorm(sys::StateSpace{Continuous, <:Dual}; kwargs...)
    A,B,C,D = ssdata(sys)
    pars = ComponentVector(; A,B,C,D)
    γ, w = implicit_hinfnorm(pars)
    γ, w
end


# ## Schur
function forward_schur(A)
    F = schur(A)
    ComponentVector(; F.Z, F.T), F
end

function conditions_schur(A, F, s)
    (; Z, T) = F
    if all(isreal, s.values)
        [
            vec(Z' * A * Z - T);
            vec(Z' * Z - I + LowerTriangular(T) - Diagonal(T))
        ]
    else
        [
            vec(Z' * A * Z - T);
            vec(Z' * Z - I + UpperTriangular(T) - Diagonal(T))
        ]
    end
end

const implicit_schur = ImplicitFunction(forward_schur, conditions_schur)

# vectors = Z
# Schur = T
# A = F.vectors * F.Schur * F.vectors'
# A = Z * T * Z'
function LinearAlgebra.schur(A::AbstractMatrix{<:Dual})
    ZT, F = implicit_schur(A)
    n = length(A)
    Z = reshape(ZT[1:n], size(A))
    T = reshape(ZT[n+1:end], size(A))
    Schur(T, Z, F.values)
end


end # module

