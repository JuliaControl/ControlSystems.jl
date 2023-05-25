using ControlSystemsBase
using ImplicitDifferentiation
using ForwardDiff
using FiniteDifferences
using ComponentArrays
fdgrad(f, x) = FiniteDifferences.grad(central_fdm(3, 1), f, x) |> first

P = ssrand(1, 1, 1)
Q = [1.0;;]
R = [1.0;;]
r = [1.0;]
q = [1.0;]


function difffun(r)
    R = reshape(r, 1, 1)
    sum(lqr(P, Q, R))
end
J1 = ForwardDiff.gradient(difffun, r)
J2 = fdgrad(difffun, r)
@test J1 ≈ J2 rtol = 1e-6

function difffun(r)
    R = reshape(r, 1, 1)
    sum(kalman(P, Q, R))
end
J1 = ForwardDiff.gradient(difffun, r)
J2 = fdgrad(difffun, r)
@test J1 ≈ J2 rtol = 1e-6


function difffun(q)
    Q = reshape(q, 1, 1)
    Rd = eltype(Q).(R)
    sum(lqr(P, Q, Rd))
end
J1 = ForwardDiff.gradient(difffun, q)
J2 = fdgrad(difffun, q)
@test J1 ≈ J2 rtol = 1e-6


P = ssrand(1, 1, 1, Ts=0.01)
function difffun(r)
    R = reshape(r, 1, 1)
    sum(lqr(P, Q, R))
end
J1 = ForwardDiff.gradient(difffun, r)
J2 = fdgrad(difffun, r)
@test J1 ≈ J2 rtol = 1e-6

function difffun(r)
    R = reshape(r, 1, 1)
    sum(kalman(P, Q, R))
end
J1 = ForwardDiff.gradient(difffun, r)
J2 = fdgrad(difffun, r)
@test J1 ≈ J2 rtol = 1e-6


function difffun(q)
    Q = reshape(q, 1, 1)
    Rd = eltype(Q).(R)
    sum(lqr(P, Q, Rd))
end
J1 = ForwardDiff.gradient(difffun, q)
J2 = fdgrad(difffun, q)
@test J1 ≈ J2 rtol = 1e-6


## Lyap
P = ssrand(1, 1, 2)
function difffun(q)
    Q = reshape(q, 2, 2)
    sum(ControlSystemsBase.lyap(P, Q))
end

q = [2.0 1; 1 2] |> vec
J1 = ForwardDiff.gradient(difffun, q)

J1 = reshape(J1, 2,2)
J1 = vec((J1 + J1') ./ 2)
J2 = fdgrad(difffun, q)
@test J1 ≈ J2 rtol = 1e-6



## hinfnorm

A = [-0.1;;]
B = [2;;]
C = [3;;]
D = 0

# A = randn(2,2) - 5I
# B = randn(2,1)
# C = randn(1,2)

sys = DemoSystems.double_mass_model(outputs=[2,4]) |> minreal
(; A, B, C, D) = sys

function difffun(pars)
    (; A,B,C,D) = pars
    hn, w = hinfnorm(ss(A, B, C, 0))
    if hn isa AbstractVector
        hn = hn[]
    end
    hn
end

pars = ComponentVector(; A,B,C,D)
J1 = ForwardDiff.gradient(difffun, pars)
J2 = fdgrad(difffun, pars)
@test J1 ≈ J2 rtol = 1e-5



## Schur currently not working when the matrix A has complex eigenvalues.

# sys = ssrand(1,1,3, proper=true)
# sys.A .-= 1I(3)
# (; A, B, C, D) = sys

# function schur_sys(sys)
#     SF = schur(sys.A)
#     A = SF.T
#     B = SF.Z'*sys.B
#     C = sys.C*SF.Z
#     ss(A,B,C,sys.D, sys.timeevol)
# end

# function difffun(pars)
#     (; A,B,C,D) = pars
#     sys = ss(A, B, C, 0)
#     sys = schur_sys(sys)
#     sum(abs2, freqresp(sys, 0.1))
# end

# pars = ComponentVector(; A,B,C,D)

# J1 = ForwardDiff.gradient(difffun, pars)[:]
# J2 = fdgrad(difffun, pars)[:]


# @show norm(J1-J2)

# @test J1 ≈ J2 rtol = 1e-3

