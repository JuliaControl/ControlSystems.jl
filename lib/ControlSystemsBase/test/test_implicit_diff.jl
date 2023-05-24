using ControlSystemsBase
using ImplicitDifferentiation
using ForwardDiff
using FiniteDifferences
using ComponentArrays

P = ssrand(1, 1, 1)
Q = [1.0;;]
R = [1.0;;]
r = [1.0;]
q = [1.0;]

gr(f, x) = FiniteDifferences.grad(central_fdm(3, 1), f, x) |> first

function difffun(r)
    R = reshape(r, 1, 1)
    sum(lqr(P, Q, R))
end
J1 = ForwardDiff.gradient(difffun, r)
J2 = gr(difffun, r)
@test J1 ≈ J2 rtol = 1e-6

function difffun(r)
    R = reshape(r, 1, 1)
    sum(kalman(P, Q, R))
end
J1 = ForwardDiff.gradient(difffun, r)
J2 = gr(difffun, r)
@test J1 ≈ J2 rtol = 1e-6


function difffun(q)
    Q = reshape(q, 1, 1)
    Rd = eltype(Q).(R)
    sum(lqr(P, Q, Rd))
end
J1 = ForwardDiff.gradient(difffun, q)
J2 = gr(difffun, q)
@test J1 ≈ J2 rtol = 1e-6


P = ssrand(1, 1, 1, Ts=0.01)
function difffun(r)
    R = reshape(r, 1, 1)
    sum(lqr(P, Q, R))
end
J1 = ForwardDiff.gradient(difffun, r)
J2 = gr(difffun, r)
@test J1 ≈ J2 rtol = 1e-6

function difffun(r)
    R = reshape(r, 1, 1)
    sum(kalman(P, Q, R))
end
J1 = ForwardDiff.gradient(difffun, r)
J2 = gr(difffun, r)
@test J1 ≈ J2 rtol = 1e-6


function difffun(q)
    Q = reshape(q, 1, 1)
    Rd = eltype(Q).(R)
    sum(lqr(P, Q, Rd))
end
J1 = ForwardDiff.gradient(difffun, q)
J2 = gr(difffun, q)
@test J1 ≈ J2 rtol = 1e-6



## hinfnorm

A = [-0.1;;]
B = [2;;]
C = [3;;]
D = 0

# A = randn(2,2) - 5I
# B = randn(2,1)
# C = randn(1,2)

sys = ssrand(1,2,4, proper=true)
(; A, B, C, D) = sys

function difffun(pars)
    (; A,B,C,D) = pars
    hn, w = hinfnorm(ss(A, B, C, 0))
    hn
end

pars = ComponentVector(; A,B,C,D)
J1 = ForwardDiff.gradient(difffun, pars)
J2 = gr(difffun, pars)
@test J1 ≈ J2 rtol = 1e-6


pars = ComponentVector(; A,B,C,D)
res, v = forward_hinfnorm(pars; tol=1e-10)

conditions_hinfnorm(pars, res, v; tol=1e-10)[]

linear_solver = (A, b) -> (Matrix(A) \ b, (solved=true,))
implicit_hinfnorm = ImplicitFunction(forward_hinfnorm, conditions_hinfnorm, linear_solver)

implicit_hinfnorm(pars)

fun = p->implicit_hinfnorm(p; tol=1e-8)[1]
g1 = ForwardDiff.jacobian(fun, pars)
g2 = gr(fun, pars)[:]'

@test g1 ≈ g2 rtol = 1e-5

