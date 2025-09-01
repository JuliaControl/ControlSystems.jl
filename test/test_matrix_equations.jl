using StaticArrays
using Test
using Printf
using Random

using MatrixEquations

Random.seed!(0)


include("../src/LyapTest.jl")

# Test _schurstructure
R = diagm(0 => [1, 1, 1, 1, 1, 1, 1, 1], -1 => [1, 0, 0, 0, 1, 0, 0])
d0 = [2, 1, 1, 2, 1, 1]
b0 = [1:2, 3:3, 4:4, 5:6, 7:7, 8:8]
@test (d0, b0, 6) == LyapTest._schurstructure(R, Val(:U))
@test (d0, b0, 6) == LyapTest._schurstructure(R', Val(:L))


# Test inner small dimension sylv solver
for T = [Float64, ComplexF64]
    for m=1:2, n=1:2
        A = randn(T, m, m)
        B = randn(T, n, n)
        C = randn(T, m, n)

        X = LyapTest._sylvd!(A, B, copy(C))

        @test Matrix(A*X*B - X - C) ≈ zeros(m,n) atol=1e-13

        X = LyapTest._sylvc!(A, B, copy(C))

        @test Matrix(A*X + X*B - C) ≈ zeros(m,n) atol=1e-13
    end
end


# Solve Sylvester equations with A = diagm(a) and B = diagm(b)
sylvc_diag = (a, b, C) -> [C[i,j]/(a[i] + b[j]) for i=1:length(a), j=1:length(b)]
sylvd_diag = (a, b, C) -> [C[i,j]/(a[i]*b[j] - 1) for i=1:length(a), j=1:length(b)]

# Compute C from X to allow convenient check of solution
sylvc_rhs = (A, B, X) -> (A*X + X*B)
sylvd_rhs = (A, B, X) -> (A*X*B - X)
lyapc_rhs = (A, X) -> -(A*X + X*A')
lyapd_rhs = (A, X) -> -(A*X*A' - X)

#
# Basic tests with diagonal 2x2 matrices
#
a = [2.0, 3]; A = diagm(a)
ac = a .+ im; Ac = diagm(ac)
b = [3.0, 5]; B = diagm(b)
C = [1.0 2; 2 1]
Cc = [1.0 2+im; 1+2im 1]
Cherm = [1.0 2+im; 2-im 1]

@test LyapTest.sylvc_schur!(A, B, copy(C)) ≈ sylvc_diag(a, b, C)
@test LyapTest.sylvc_schur!(Ac, B, copy(Cc)) ≈ sylvc_diag(ac, b, Cc)
@test LyapTest.sylvc(A, B, C) ≈ sylvc_diag(a, b, C)
@test LyapTest.sylvc(Ac, B, Cc) ≈ sylvc_diag(ac, b, Cc)

@test LyapTest.sylvd_schur!(A, B, copy(C)) ≈ sylvd_diag(a, b, C)
@test LyapTest.sylvd_schur!(Ac, B, copy(Cc)) ≈ sylvd_diag(ac, b, Cc)
@test LyapTest.sylvd(A, B, C) ≈ sylvd_diag(a, b, C)
@test LyapTest.sylvd(Ac, B, Cc) ≈ sylvd_diag(ac, b, Cc)

@test LyapTest.lyapc_schur!(A, copy(C)) ≈ sylvc_diag(a, a, -C)
@test LyapTest.lyapc_schur!(Ac, copy(Cherm)) ≈ sylvc_diag(conj(ac), ac, -Cherm)
@test LyapTest.lyapc(A, C) ≈ sylvc_diag(a, a, -C)
@test LyapTest.lyapc(Ac, Cherm) ≈ sylvc_diag(conj(ac), ac, -Cherm)

@test LyapTest.lyapd_schur!(A, copy(C)) ≈ sylvd_diag(a, a, -C)
@test LyapTest.lyapd_schur!(Ac, copy(Cherm)) ≈ sylvd_diag(ac, conj(ac), -Cherm)
@test LyapTest.lyapd(A, C) ≈ sylvd_diag(a, a, -C)
@test LyapTest.lyapd(Ac, Cherm) ≈ sylvd_diag(ac, conj(ac), -Cherm)

# A few tests for the naive version, should have more
@test LyapTest.sylvc(A, B, C, Val(:naive)) ≈ sylvc_diag(a, a, -C)
@test LyapTest.lyapc(A, C, Val(:naive)) ≈ sylvc_diag(a, a, -C)



#
# Tests with small trinagular matrices of various types
#
A3 = [2 0 0; 3 4 0; 5 6 7]
B3 = [2 3 4; 0 1 1; 0 0 2]
X3 = reshape(1:9, 3, 3)

for (A, B, X, tsname) in [(fill(2.0, 1, 1), fill(3.0, 1, 1), fill(1.0, 1, 1), "1x1"),
                  (diagm([2.0,3]), diagm([4.0,5]), [1.0 2; 3 4], "2x2 diagonal"),
                  ([-2 0; 2 3], [3.0 4; 0 5], [1.0 2; 3 4], "2x2 tringular"),
                  ([-2.0 0; 2 3], fill(5.0, 1, 1), [1.0; 2][:,:], "2x1"),
                  ([-2 0; 2 3], fill(5, 1, 1), [1; 2][:,:], "2x1"),
                  (big.([-2.0 0; 2 3]), fill(5//2, 1, 1), [1; 2][:,:], "2x1"),
                  (fill(2.0, 1, 1), [3.0 4; 0 5], [1.0 2], "1x2"),
                  (big.(A3), Rational.(B3), X3, "3x3"),
                  ]

    # Generate complex version, while keeping the structure
    Ac = A .* (1 + im)
    Bc = B .* (1 + im)
    Xc = X .* (1 + im)

    @testset "sylv(c/d)_schur!, $(tsname) $(eltype.((A, B, X)))" begin
        @test LyapTest.sylvc_schur!(A, B, sylvc_rhs(A, B, X)) ≈ X
        @test LyapTest.sylvc_schur!(Ac, B, sylvc_rhs(Ac, B, X)) ≈ X
        @test LyapTest.sylvc_schur!(A, B, sylvc_rhs(A, B, Xc)) ≈ Xc
        @test LyapTest.sylvc_schur!(Ac, Bc, sylvc_rhs(Ac, Bc, Xc)) ≈ Xc

        @test LyapTest.sylvd_schur!(A, B, sylvd_rhs(A, B, X)) ≈ X
        @test LyapTest.sylvd_schur!(Ac, B, sylvd_rhs(Ac, B, X)) ≈ X
        @test LyapTest.sylvd_schur!(A, B, sylvd_rhs(A, B, Xc)) ≈ Xc
        @test LyapTest.sylvd_schur!(Ac, Bc, sylvd_rhs(Ac, Bc, Xc)) ≈ Xc

    end

    if size(X,1) == size(X,2)
        Xherm = X + X'
        Xcherm = Xc + Xc'

        @testset "lyap(c/d)_schur!, $(tsname) $(eltype.((A, C)))" begin
            @test LyapTest.lyapc_schur!(A, lyapc_rhs(A, Xherm)) ≈ Xherm
            @test LyapTest.lyapc_schur!(Ac, lyapc_rhs(Ac, Xherm)) ≈ Xherm
            @test LyapTest.lyapc_schur!(A, lyapc_rhs(A, Xcherm)) ≈ Xcherm
            @test LyapTest.lyapc_schur!(Ac, lyapc_rhs(Ac, Xcherm)) ≈ Xcherm

            @test LyapTest.lyapd_schur!(A, lyapd_rhs(A, Xherm)) ≈ Xherm
            @test LyapTest.lyapd_schur!(Ac, lyapd_rhs(Ac, Xherm)) ≈ Xherm
            @test LyapTest.lyapd_schur!(A, lyapd_rhs(A, Xcherm)) ≈ Xcherm
            @test LyapTest.lyapd_schur!(Ac, lyapd_rhs(Ac, Xcherm)) ≈ Xcherm
        end
    end
end


Random.seed!(0)
for (m,n) in [(1,1), (1,2), (2,1), (2,2), (3,3), (3,2), (2,3), (5,5), (7,8), (100,100)]
    for (TA, TB, TC) in [(Float64, Float64, Float64),
                         (ComplexF64, Float64, Float64),
                         (Float64, ComplexF64, Float64),
                         (Float64, ComplexF64, ComplexF64),
                         (ComplexF64, ComplexF64, ComplexF64),
                         (Int64, ComplexF64, ComplexF64),
                         (Int64, Int64, Int64)]

        # Generate complex version, while keeping the structure
        A = TA <: Int ? rand(1:50, m, m) : randn(TA, m, m)
        B = TB <: Int ? rand(1:50, n, n) : randn(TB, n, n)
        C = TC <: Int ? rand(1:50, m, n) : randn(TC, m, n)

        rtol = m > 50 ? 1e-10 : 1e-12

        @testset "lyap/sylv ($m, $n, $TA, $TB, $TC)" begin
            X =  LyapTest.sylvc(A, B, C)
            @test A*X + X*B ≈ C rtol=rtol

            X =  LyapTest.sylvd(A, B, C)
            @test A*X*B - X ≈ C rtol=rtol

            if m == n
                Q = Matrix(Hermitian(randn(m,m)))

                X = LyapTest.lyapc(A, Q)
                @test A*X + X*A' ≈ -Q rtol=rtol

                X = LyapTest.lyapd(A, Q)
                @test A*X*A' - X ≈ -Q rtol=rtol
            end
        end
    end
end
