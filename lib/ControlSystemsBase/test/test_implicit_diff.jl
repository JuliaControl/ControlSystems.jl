using ControlSystemsBase
using ImplicitDifferentiation
using ForwardDiff
using FiniteDifferences
using ComponentArrays
using Test, LinearAlgebra
fdgrad(f, x) = FiniteDifferences.grad(central_fdm(3, 1), f, x) |> first

P = ssrand(1, 1, 2)
Q = [1.0 0.1; 0.1 0.5]
R = [1.0;;]
r = [1.0;]
q = vec(Q)


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
    Q = reshape(q, 2, 2)
    Q = (Q .+ Q') ./ 2 # Needed for finite diff
    Rd = eltype(Q).(R)
    sum(lqr(P, Q, Rd))
end
J1 = ForwardDiff.gradient(difffun, q)
J2 = fdgrad(difffun, q)
@test J1 ≈ J2 rtol = 1e-6


P = ssrand(1, 1, 2, Ts=0.01)
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
    Q = reshape(q, 2, 2)
    Q = (Q .+ Q') ./ 2 # Needed for finite diff
    Rd = eltype(Q).(R)
    sum(lqr(P, Q, Rd))
end
J1 = ForwardDiff.gradient(difffun, q)
J2 = fdgrad(difffun, q)
@test J1 ≈ J2 rtol = 1e-6


## Keyword arguments and output sizes of kalman and lqr (ControlSystems.jl issue 1105)
dualize(R) = ForwardDiff.Dual{:test}.(R, 1.0)
let
    for Ts in (nothing, 0.01)
        P = Ts === nothing ? ssrand(2, 1, 3) : ssrand(2, 1, 3; Ts)
        Q = [2.0 0.1 0; 0.1 2 0; 0 0 1]
        R = [1.0 0.2; 0.2 1.5]
        r = vec(R)
        directs = Ts === nothing ? (false,) : (false, true)
        for direct in directs
            Kref, Xref, pref = kalman(P, Q, R; direct, extra=Val(true))
            Kd, Xd, pd = kalman(P, Q, dualize(R); direct, extra=Val(true))
            @test size(Kd) == size(Kref) == (P.nx, P.ny)
            @test ForwardDiff.value.(Kd) ≈ Kref
            @test ForwardDiff.value.(Xd) ≈ Xref
            @test sort(pd, by=x->(real(x), imag(x))) ≈ sort(pref, by=x->(real(x), imag(x)))
            @test size(kalman(P, Q, dualize(R); direct)) == (P.nx, P.ny)

            difffun = function (r)
                Rd = reshape(r, 2, 2)
                Rd = (Rd .+ Rd') ./ 2
                K = kalman(P, Q, Rd; direct)
                sum(abs2, K * [1.0 2.0; 3.0 4.0])
            end
            @test ForwardDiff.gradient(difffun, r) ≈ fdgrad(difffun, r) rtol = 1e-5

            difffun_extra = function (r)
                Rd = reshape(r, 2, 2)
                Rd = (Rd .+ Rd') ./ 2
                K, X, _ = kalman(P, Q, Rd; direct, extra=Val(true))
                sum(abs2, K * [1.0 2.0; 3.0 4.0]) + sum(X)
            end
            @test ForwardDiff.gradient(difffun_extra, r) ≈ fdgrad(difffun_extra, r) rtol = 1e-5
        end
        Ts === nothing && @test_throws ErrorException kalman(P, Q, dualize(R); direct=true)

        Pu = Ts === nothing ? ssrand(1, 2, 3) : ssrand(1, 2, 3; Ts)
        Kref, Sref, pref = lqr(Pu, Q, R; extra=Val(true))
        Kd, Sd, pd = lqr(Pu, Q, dualize(R); extra=Val(true))
        @test size(Kd) == size(Kref) == (Pu.nu, Pu.nx)
        @test ForwardDiff.value.(Kd) ≈ Kref
        @test ForwardDiff.value.(Sd) ≈ Sref
        @test sort(pd, by=x->(real(x), imag(x))) ≈ sort(pref, by=x->(real(x), imag(x)))
    end

    # Example from the issue, with coupled state variables to obtain a nonzero derivative
    function cost(p::T) where T <: Real
        A = [0.9 0.1; 0.0 p]
        C = [1.0 0.0]
        R1 = Matrix{T}(I, 2, 2)
        R2 = Matrix{T}(I, 1, 1)
        K = kalman(Discrete, A, C, R1, R2, direct=true)
        @test size(K) == (2, 1)
        K2, _ = kalman(Discrete, A, C, R1, R2, extra=Val(true))
        @test size(K2) == (2, 1)
        sum(K) + sum(K2)
    end
    @test ForwardDiff.derivative(cost, 0.1) ≈ central_fdm(5, 1)(cost, 0.1) rtol = 1e-6
end


# Diff w.r.t. plant
P = ssrand(1, 1, 2)
Q = I(P.nx)
R = [1.0;;]
r = [1.0;]

function difffun(a)
    A = reshape(a, 2, 2)
    Rd = eltype(a).(R)
    sum(lqr(Continuous, A, P.B, Q, Rd))
end
a = P.A[:]
J1 = ForwardDiff.gradient(difffun, a)
J2 = fdgrad(difffun, a)
@test J1 ≈ J2 rtol = 1e-6

## are directly
Q = [2.0 0.1; 0.1 2]
q = vec(Q)
function difffun(q)
    Q = reshape(q, 2, 2)
    Q = (Q .+ Q') ./ 2 # Needed for finite diff
    Rd = eltype(Q).(R)
    sum(are(P, Q, Rd))
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

## positive definite Lyap
Ql = [1 0; 0.1 1]
function difffun(q)
    Ql = LowerTriangular(copy(reshape(q, 2, 2)))
    sum(ControlSystemsBase.plyap(P, Ql))
end

q = Ql |> vec
J1 = ForwardDiff.gradient(difffun, q)
J2 = fdgrad(difffun, q)
@test J1 ≈ J2 rtol = 1e-6



## covar (tests plyap)
P = ssrand(1, 2, 2, proper=true)
function difffun(q)
    Q = reshape(q, 2, 2)
    Q = (Q .+ Q') ./ 2 # Needed for finite diff
    sum(ControlSystemsBase.covar(P, Q))
end

q = Q |> vec
J1 = ForwardDiff.gradient(difffun, q)
J2 = fdgrad(difffun, q)
@test J1 ≈ J2 rtol = 1e-6

# covar w.r.t. plant
function difffun(a)
    A = copy(reshape(a, 2, 2))
    sys2 = ss(A, P.B, P.C, P.D)
    sum(ControlSystemsBase.covar(sys2, Q))
end

a = P.A |> vec
J1 = ForwardDiff.gradient(difffun, a)
J2 = fdgrad(difffun, a)
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
    if hn isa AbstractArray
        hn = hn[]
    end
    hn
end

pars = ComponentVector(; A,B,C,D)
J1 = ForwardDiff.gradient(difffun, pars)
J2 = fdgrad(difffun, pars)
@test J1 ≈ J2 rtol = 1e-5



## Schur decomposition

sys = ssrand(1,1,3, proper=true)
# sys.A .-= 1I(3)
(; A, B, C, D) = sys

function schur_sys(sys)
    SF = schur(sys.A)
    A = SF.T
    B = SF.Z'*sys.B
    C = sys.C*SF.Z
    ss(A,B,C,sys.D, sys.timeevol)
end

function difffun(pars)
    (; A,B,C,D) = pars
    sys = ss(A, B, C, 0)
    sys = schur_sys(sys)
    sum(abs2, freqresp(sys, 0.1))
end

pars = ComponentVector(; A,B,C,D)

J1 = ForwardDiff.gradient(difffun, pars)[:]
J2 = fdgrad(difffun, pars)[:]
# @show norm(J1-J2)
@test J1 ≈ J2 rtol = 1e-5

## Test differentiation of tf 2 ss conversion

function difffun(pars)
    P = tf(1, pars)
    sum(abs2, step(P, 0:0.1:10, method=:tustin).y + 
        impulse(P, 0:0.1:10, method=:tustin).y
    )
end

pars = [1.0, 2, 1]

g1 = ForwardDiff.gradient(difffun, pars)
g2 = fdgrad(difffun, pars)
@test g1 ≈ g2 rtol = 1e-5