module ControlSystemsBaseImplicitDifferentiationExt
using ControlSystemsBase
import ControlSystemsBase: lqr, are, ContinuousType, DiscreteType
using ForwardDiff
using ForwardDiff: Dual
using ImplicitDifferentiation
using ComponentArrays
using LinearAlgebra


# The solver keyword arguments `kw` are a constant argument of the implicit function. The
# byproduct contains the remaining outputs of the Riccati solver (closed-loop eigenvalues, gain etc.),
# which are not differentiated.
function forward_arec(pars, kw)
    (; A,B,Q,R) = pars
    X, rest... = ControlSystemsBase.arec(A, B, R, Q; kw...)
    X, rest
end

function conditions_arec(pars, X, rest, kw)
    (; A,B,Q,R) = pars
    C = A'X
    C .+= C'
    C .+= Q
    XB = X*B
    mul!(C, XB, R\XB', -1, 1)
end

const implicit_arec = ImplicitFunction(forward_arec, conditions_arec)

function forward_ared(pars, kw)
    (; A,B,Q,R) = pars
    X, rest... = ControlSystemsBase.ared(A, B, R, Q; kw...)
    X, rest
end

function conditions_ared(pars, X, rest, kw)
    # A'XA - X - (A'XB)(R+B'XB)^(-1)(B'XA) + Q = 0
    (; A,B,Q,R) = pars
    AX = A'X
    C = AX*A
    C .+= Q .- X'
    AXB = AX*B
    C .-= AXB*(((B'X*B) .+= R)\AXB')
    C
end

const implicit_ared = ImplicitFunction(forward_ared, conditions_ared)

# Returns the solution X, differentiable with respect to all of A, B, Q, R, and the remaining
# outputs of the Riccati solver, starting with the closed-loop eigenvalues, which carry no
# derivative information.
function _are_extra(te, A, B, Q, R; kwargs...)
    pars = ComponentVector(; A,B,Q,R)
    implicit = te isa ContinuousType ? implicit_arec : implicit_ared
    X0, rest = implicit(pars, NamedTuple(kwargs))
    reshape(X0, size(A)), rest
end

"""
    are(::Continuous, A, B, Q, R::AbstractMatrix{<:Dual}; kwargs...)
    are(::Discrete,   A, B, Q, R::AbstractMatrix{<:Dual}; kwargs...)

To make the ARE solver work with dual numbers, make sure that the `R` matrix has the dual number element type. The keyword arguments are passed to the Riccati solver.
"""
function are(te::ContinuousType, A::AbstractMatrix, B, Q, R::AbstractMatrix{<:Dual}; kwargs...)
    first(_are_extra(te, A, B, Q, R; kwargs...))
end

function are(te::DiscreteType, A::AbstractMatrix, B, Q, R::AbstractMatrix{<:Dual}; kwargs...)
    first(_are_extra(te, A, B, Q, R; kwargs...))
end

"""
    lqr(::Continuous, A, B, Q, R::AbstractMatrix{<:Dual}; extra = Val(false), kwargs...)
    lqr(::Discrete,   A, B, Q, R::AbstractMatrix{<:Dual}; extra = Val(false), kwargs...)

To make the LQR solver work with dual numbers, make sure that the `R` matrix has the dual number element type. This method is also used by [`kalman`](@ref) when `R2` has the dual number element type.

If `extra = Val(true)`, the outputs are the same as for the generic method. The closed-loop eigenvalues `p` and the additional outputs of the Riccati solver carry no derivative information.
"""
function lqr(te::ContinuousType, A, B, Q, R::AbstractMatrix{<:Dual}; extra::Val{E} = Val(false), kwargs...) where E
    X, (p, _, rest...) = _are_extra(te, A, B, Q, R; kwargs...)
    K = R\(B'X)
    E ? (K, X, p, rest...) : K
end

function lqr(te::DiscreteType, A, B, Q, R::AbstractMatrix{<:Dual}; extra::Val{E} = Val(false), kwargs...) where E
    X, (p, _, rest...) = _are_extra(te, A, B, Q, R; kwargs...)
    BX = B'X
    K = (R+BX*B)\(BX*A)
    E ? (K, X, p, rest...) : K
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
    O .+ (X .- X')
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
    Matrix(ControlSystemsBase.plyapc(A, Q)), 0
end

function conditions_plyapc(pars, Xc, noneed)
    (; A,Q) = pars
    Q = Q*Q'
    X = Xc*Xc'
    AX = A*X
    O = AX .+ AX' .+ Q
    O .+ (Xc .- UpperTriangular(Xc))
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
    γ, w = hinfnorm(sys; kwargs...)
    return [γ], w
end

function conditions_hinfnorm(pars, γ_vec, w; tol=1e-10)
    γ = only(γ_vec)
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
    γ_vec, w = implicit_hinfnorm(pars)
    only(γ_vec), w
end


# ## Schur
function forward_schur(A)
    F = schur(A)
    ComponentVector(; F.Z, F.T), F
end

function conditions_schur(A, F, s)
    (; Z, T) = F
    if all(isreal, s.values)
        ComponentVector(; Z = Z' * A * Z - T, T = Z' * Z - I + LowerTriangular(T) - Diagonal(T))
    else
        ComponentVector(; Z = Z' * A * Z - T, T = Z' * Z - I + UpperTriangular(T) - Diagonal(T))
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
