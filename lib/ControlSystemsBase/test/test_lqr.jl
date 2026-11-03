using StaticArrays: SMatrix, SVector

@testset "LQR cross term" begin
    @testset "scalar reference solutions" begin
        A = B = N = ones(1, 1)
        K, P = lqr(Continuous, A, B, 2I, 2I, N; extra=Val(true))
        @test K ≈ fill(2.0, 1, 1)
        @test P ≈ fill(3.0, 1, 1)
        K, P = lqr(Discrete, A, B, 2I, 2I, N; extra=Val(true))
        @test K ≈ fill(sqrt(3) - 1, 1, 1)
        @test P ≈ fill(sqrt(3), 1, 1)
    end

    @testset "$te, $T" for te in (Continuous, Discrete), T in (Float32, Float64, ComplexF64)
        A = T[0 1 0; 0 0 1; -1 -2 -3]
        B = T[0 1; 1 0; 1 1]
        N = T[0.2 -0.1; 0.1 0.3; -0.2 0.1]
        R = T[2 0.3; 0.3 1]
        if T <: Complex
            A[1, 2] += 0.2im
            B[2, 1] += 0.1im
            N[1, 2] += 0.1im
            R[1, 2] += 0.1im
            R[2, 1] -= 0.1im
        end
        Q = Diagonal(T[2, 3, 4])
        original = deepcopy((A, B, Q, R, N))
        result = lqr(te, A, B, Q, R, N; extra=Val(true))
        K, P, p = result
        tol = 100eps(real(T))

        @test size(K) == (2, 3)
        @test eltype(K) == T
        @test P ≈ P'
        @test lqr(te, A, B, Q, R, N) ≈ K
        @test approxsetequal(p, eigvals(A - B*K); atol=tol, rtol=tol)
        if te == Continuous
            @test all(real.(p) .< 0)
            @test K ≈ R \ (B'*P + N')
            scale = 2norm(A)*norm(P) + norm(P*B + N)*norm(K) + norm(Q)
            @test norm(A'*P + P*A - (P*B + N)*K + Q) <= tol*scale
        else
            @test all(abs.(p) .< 1)
            @test K ≈ (R + B'*P*B) \ (B'*P*A + N')
            scale = norm(A)^2*norm(P) + norm(P) + norm(A'*P*B + N)*norm(K) + norm(Q)
            @test norm(A'*P*A - P - (A'*P*B + N)*K + Q) <= tol*scale
        end

        # Complete the square using v = u + F*x, in either time domain.
        F = R \ N'
        Abar = A - B*F
        Qbar = Hermitian(Q - N*F)
        Kbar, Pbar, pbar = lqr(te, Abar, B, Qbar, R; extra=Val(true))
        @test K ≈ Kbar + F
        @test P ≈ Pbar
        @test approxsetequal(p, pbar; atol=tol, rtol=tol)

        sys = te == Continuous ? ss(A, B, I, 0) : ss(A, B, I, 0, 0.1)
        @test lqr(sys, Q, R, N) ≈ K
        @test lqr(sys, Q, R, N; extra=Val(true))[2] ≈ P
        @test lqr(sys.timeevol, A, B, Q, R, N) ≈ K
        @test lqr(te, A, B, Q, R, N; rtol=eps(real(T))) ≈ K
        @test lqr(te, A, B, Q, R, zeros(T, size(N))) ≈ lqr(te, A, B, Q, R)
        @test lqr(te, A, B, Q, R, -N) ≉ K
        @test (A, B, Q, R, N) == original

        # Preserve all backend outputs when extra=Val(true), including solver metadata.
        solver = te == Continuous ? ControlSystemsBase.MatrixEquations.arec : ControlSystemsBase.MatrixEquations.ared
        expected = solver(A, B, R, Q, N)
        @test length(result) == length(expected)
        @test all(isequal.(result[4:end], expected[4:end]))
    end

    @testset "single input and matrix representations ($te)" for te in (Continuous, Discrete)
        A = [0.0 1; 0 0] # Also exercises a singular discrete-time state matrix.
        B = [0.0; 1]
        N = [0.2; -0.1]
        K = lqr(te, A, B, I, I, N)
        @test lqr(te, A, reshape(B, 2, 1), I, I, reshape(N, 2, 1)) ≈ K
        @test lqr(te, A, B, 1, 1, N) ≈ K
        @test lqr(te, A, B, Diagonal(ones(2)), Diagonal(ones(1)), N) ≈ K
        @test lqr(te, A, B, Symmetric(Matrix(1.0I, 2, 2)), Hermitian(ones(1, 1)), N) ≈ K
        @test lqr(te, view(A, :, :), view(B, :), I, I, view(N, :)) ≈ K
        @test lqr(te, A, B, I, I, reshape(N, 1, 2)') ≈ K
        @test lqr(te, A, B, I, I, SVector{2}(N)) ≈ K
        @test lqr(te, A, reshape(B, 2, 1), I, I, SMatrix{2,1}(N)) ≈ K
        @test lqr(te, A, B, I, I, sparse(reshape(N, 2, 1))) ≈ K
    end

    @testset "semidefinite transformed cost ($te)" for te in (Continuous, Discrete)
        B = Matrix(1.0I, 2, 2)
        N = Diagonal([0.2, -0.1])
        R = Diagonal([2.0, 1.0])
        F = R \ N'
        Q = N*F
        A = (te == Continuous ? -B : 0.5B) + B*F
        K, P, p = lqr(te, A, B, Q, R, N; extra=Val(true))
        @test K ≈ F
        @test norm(P) < 1e-14
        @test approxsetequal(p, eigvals(A - B*F); atol=1e-14)
    end

    @testset "dimensions ($te)" for te in (Continuous, Discrete)
        A = [0.0 1 0; 0 0 1; -1 -2 -3]
        B = [0.0 1; 1 0; 1 1]
        for N in (zeros(2, 2), zeros(3, 1), zeros(2, 3), zeros(3), zeros(3, 3))
            @test_throws DimensionMismatch lqr(te, A, B, I, I, N)
        end
        @test_throws DimensionMismatch lqr(te, A, B[1:2, :], I, I, zeros(3, 2))
        @test_throws DimensionMismatch lqr(te, A, B, ones(2, 2), I, zeros(3, 2))
        @test_throws DimensionMismatch lqr(te, A, B, I, ones(3, 3), zeros(3, 2))
        @test_throws DimensionMismatch lqr(te, A, B[:, 1], I, I, zeros(3, 2))
    end
end
