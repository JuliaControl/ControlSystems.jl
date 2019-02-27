using ControlSystems
using Test
using LinearAlgebra
using Random

@testset "H-infinity design" begin
  """
  Tests for the public and private methods of the hInfSynthesis function. This
  function utilizes the preexisting ControlSystems toolbox, and performs a
  H-infinity synthesis using the dual Riccati equation approach. As such,
  the synthesis is done in a set of steps.

  (1) Re-writing the specifications on an extended state-space form.
  (2) Verifying that the resulting extended state-space object satisfies a set of
      assumptions required for proceeding with the synthesis.
  (3) A coordinate transform to enable the synthesis.
  (4) Synthesis using the γ-iterations, checking if a solution to the H-infinity
      problem exists in each iteration and applying a bisection method.
  (5) Re-transforming the system to the original coordinates
  (6) Verification that the computed solution is correct.

  In addition to these six ponts, the code also enables

  (7) A bilinear discretization with an inverse operation to move from continuous
  to discrete time, thereby enabling approximate discrete-time synthesis.
  (8) Plotting functionality to visualize the H-infinity synthesis.
  (9) Three examples which can be used to demonstrate the tool.

  Many of the tests are quite intuitive, and all points (1)-(9) are tested
  extensively with detailed comments for each test-set.
  """

  @testset "(2) Assumptions" begin
    """
    Tests the methods used to check that the assumptions are satisfied,
    incorporating two separate algorithms to determine detectability and
    observability, and a rectangular pseudoinverse to check the rank conditions.
    """

    @testset "Stabilizability check" begin
      """
      Test the check for stabilizability using the Hautus Lemma
      """
      ### Fixture
      Random.seed!(0); N = 10; M = 5;
      R = rand(Float64, (N,N))
      Q = eigvecs(R + R');
      L = rand(Float64,N).-0.5;
      A = Q*Diagonal(L)*Q';
      E = eigvals(A);

      # Here we should return false if for any positive eigenvalue λ of A, we get
      # that rank(A-λI, B) < N. We define B as conlums which are linearly
      # dependent with columns in A-λI, and check that the system returns false
      # in every case where λ > 0 and never when λ≦0.
      for ii = 1:N
        for jj = 1:(N-M)
          Ahat = A - E[ii]*Matrix{Float64}(I, N, N);
          B    = Ahat[:,jj:(jj+M)]
          if E[ii] >= 0
            @test !ControlSystems._is_stabilizable(A,B)
          else
            @test ControlSystems._is_stabilizable(A,B)
          end
        end
      end

      # Check common input types and ensure that a method error is thrown when
      # not using two abstract matrices, but some other common type
      N = 10; M = 5;
      A = rand(Float64, (N,N));
      B = rand(Float64, (N,M));
      @test_throws MethodError ControlSystems._is_detectable(A,nothing)
      @test_throws MethodError ControlSystems._is_detectable(nothing,B)
      @test_throws MethodError ControlSystems._is_detectable(A,[])
      @test_throws MethodError ControlSystems._is_detectable([],B)
      @test_throws MethodError ControlSystems._is_detectable(A,ss(1))
      @test_throws MethodError ControlSystems._is_detectable(ss(1),B)
      @test_throws MethodError ControlSystems._is_detectable(A,tf(1))
      @test_throws MethodError ControlSystems._is_detectable(tf(1),B)
    end

    @testset "Detectability check" begin
      """
      Test the check for detectability using the Hautus Lemma
      """
      ### Fixture
      Random.seed!(0); N = 10; M = 5;
      R = rand(Float64, (N,N))
      Q = eigvecs(R + R');
      L = rand(Float64,N).-0.5;
      A = Q*Diagonal(L)*Q';
      E = eigvals(A);

      # Here we should return false if for any positive eigenvalue λ of A, we get
      # that rank(A-λI; C) < N. We define C as rows which are linearly
      # dependent with rows in A-λI, and check that the system returns false
      # in every case where λ > 0 and never when λ≦0.
      for ii = 1:N
        for jj = 1:(N-M)
          Ahat = A - E[ii]*Matrix{Float64}(I, N, N);
          C    = Ahat[jj:(jj+M),:]
          if E[ii] >= 0
            @test !ControlSystems._is_detectable(A,C)
          else
            @test ControlSystems._is_detectable(A,C)
          end
        end
      end

      # Check common input types and ensure that a method error is thrown when
      # not using two abstract matrices, but some other common type
      N = 10; M = 5;
      A = rand(Float64, (M,M));
      C = rand(Float64, (N,M));
      @test_throws MethodError ControlSystems._is_detectable(A,nothing)
      @test_throws MethodError ControlSystems._is_detectable(nothing,C)
      @test_throws MethodError ControlSystems._is_detectable(A,[])
      @test_throws MethodError ControlSystems._is_detectable([],C)
      @test_throws MethodError ControlSystems._is_detectable(A,ss(1))
      @test_throws MethodError ControlSystems._is_detectable(ss(1),C)
      @test_throws MethodError ControlSystems._is_detectable(A,tf(1))
      @test_throws MethodError ControlSystems._is_detectable(tf(1),C)
    end

    @testset "Pseudoinverse computation" begin
      """
      The g-inverse is used to check the rank conditions across all frequencies,
      see e.g. assumptions A5 and A6 in the reference [1]. As was shown in the
      complementary report [2], this can be done by computing a pseudoinverse
      and checking its rank. This code tests that the compuatation of these
      g-inverses are done correctly for arbitrary matrices
      """

      ### Fixture
      tolerance = 1e-10;
      R = rand(Float64, (10,10))
      Q = eigvecs(R + R');

      # Check that the g-invere works for a set of full-rank square and
      # rectangular matrices
      for M = 1:5:11
        for N = 1:5:11
          AMN = rand(M,N);
          Pinv, status = ControlSystems._compute_pseudoinverse(AMN);
          @test status
          if M < N
            @test opnorm(AMN*Pinv - Matrix{Float64}(I, min(M,N), min(M,N))) < tolerance
          else
            @test opnorm(Pinv*AMN - Matrix{Float64}(I, min(M,N), min(M,N))) < tolerance
          end
        end
      end

      # Check common input types and ensure that a method error is thrown when
      # not using two abstract matrices, but some other common type
      @test_throws MethodError ControlSystems._compute_pseudoinverse(nothing)
      @test_throws MethodError ControlSystems._compute_pseudoinverse([])
      @test_throws MethodError ControlSystems._compute_pseudoinverse(ss(1))
      @test_throws MethodError ControlSystems._compute_pseudoinverse(tf(1))

      # Check that the method correctly reports that no pseudoinverse exists
      # when the matrix M is rank deficient
      for M = 1:5:11
        for N = 1:5:11
          U,S,V = svd(rand(M,N))
          S[1] = 0.0;
          AMN_rank_deficient = U*Diagonal(S)*V';
          Pinv, status = ControlSystems._compute_pseudoinverse(AMN_rank_deficient);
          @test !status
          @test isa(Pinv, Array{Any,1})
        end
      end
    end

    # TODO: write tests using the above submethods directly in hInf_assumptions
  end

  @testset "(4) Gamma iterations" begin
    """
    Tests the core methods of the gamma-iteration bisection method, including
    the ARE hamiltonian solver, the eigenvalue solver, the feasibility checks.
    """

    @testset "Solution feasibility check" begin
      """
      Check that a solution to the dual riccati equations is correctly reported
      as being feasible if satisfying the conditions of positive definiteness on
      the X and Y solutions are met, and the spectral radius condition of X*Y is
      also met
      """
      # Fixture
      Random.seed!(0);
      tolerance = 1e-10;
      iteration = 1;
      Random.seed!(0); N = 10; M = 5;
      R = rand(Float64, (N,N))
      Q = eigvecs(R + R');
      ρX = rand()
      ρY = rand()
      LX= rand(Float64,N); LX = ρX*sort(LX / maximum(LX))
      LY= rand(Float64,N); LY = ρY*sort(LY / maximum(LY))
      Xinf = Q*Diagonal(LX)*Q';
      Yinf = Q*Diagonal(LY)*Q';
      gamma = 1;

      # Test that the fesibility is true if ρ(Xinf*Yinf) < γ^2, that is the
      # check should be true for any
      #
      #    γ = sqrt(ρX*ρY) + ϵ
      #
      # with any epsilon greater than or equal to zero.
      @test !ControlSystems._checkFeasibility(Xinf, Yinf, sqrt(ρX*ρY)-tolerance, tolerance, iteration; verbose=false)
      @test !ControlSystems._checkFeasibility(Xinf, Yinf, sqrt(ρX*ρY),           tolerance, iteration; verbose=false)
      @test  ControlSystems._checkFeasibility(Xinf, Yinf, sqrt(ρX*ρY)+tolerance, tolerance, iteration; verbose=false)

      # Test that errors are thrown if the matrix Xinf and Yinf are not PSD down
      # to the numerical tolerance.
      L = LX;
      L[1] += -L[1] + 2*tolerance; Xpos   = Q*Diagonal(L)*Q'; # slightly positive eigenvalue
      L[1] += -L[1]; Xzero  = Q*Diagonal(LX)*Q';               # exactly one zero eigenvalue
      L[1] += -L[1] - 2*tolerance; Xneg = Q*Diagonal(L)*Q';   # slightly negative eigenvalue
      @test  ControlSystems._checkFeasibility(Xpos,  Yinf, sqrt(ρX*ρY)+tolerance, tolerance, iteration; verbose=false)
      @test  ControlSystems._checkFeasibility(Xzero, Yinf, sqrt(ρX*ρY)+tolerance, tolerance, iteration; verbose=false)
      @test !ControlSystems._checkFeasibility(Xneg,  Yinf, sqrt(ρX*ρY)+tolerance, tolerance, iteration; verbose=false)

      L = LY;
      L[1] += -L[1] + 2*tolerance; Ypos   = Q*Diagonal(L)*Q'; # slightly positive eigenvalue
      L[1] += -L[1]; Yzero  = Q*Diagonal(L)*Q';               # exactly one zero eigenvalue
      L[1] += -L[1] - 2*tolerance; Yneg = Q*Diagonal(L)*Q';   # slightly negative eigenvalue
      @test  ControlSystems._checkFeasibility(Xinf, Ypos,  sqrt(ρX*ρY)+tolerance, tolerance, iteration; verbose=false)
      @test  ControlSystems._checkFeasibility(Xinf, Yzero, sqrt(ρX*ρY)+tolerance, tolerance, iteration; verbose=false)
      @test !ControlSystems._checkFeasibility(Xinf, Yneg,  sqrt(ρX*ρY)+tolerance, tolerance, iteration; verbose=false)
    end

    # TODO: Include a check to verify that the bisection works as intended.
    # TODO: Check to verify that the Hamiltonian Shur-solver is working
    # TODO: Check to verify that the Hamiltonian Eigenvalue-solver is working
  end

  @testset "Bilinear discretization" begin
    """
    This tests the bilinear method of discretizing a continuous time StateSpace
    or ExtendedStateSpace object, moving from the Laplace-domain to the Z-domain.
    However, importantly, the method also enables approximating discrete-time
    systems with continuous-time state-space domain, moving from the Z-domain back
    to the laplace-domain. Since the L∞-norm is invariant over the bilinear
    discretization and "continuization", this effectively allows us to approximate
    discrete-time systems with a continuous-time equivalent, design an H∞-optimal
    controller for the continuous-time plant, and then re-discretizing this
    controller.

    The structure of these tests is to simply discretize a set of plants with
    various dimensions, transport them back into continuous-time domain, and then
    to discrete-time again. By then comparing the resulting two pairs of discrete-
    time and continuous-time systems in the frequency domain, we can test if the
    bilinear discretization operates as expected.
    """
    @testset "SS data" begin
      # Fixture for the SS-type data
      Random.seed!(0);
      N = [1,5,10]
      M = [1,5,10]
      H = [0.1,0.01,0.001]
      tolerance = 1e-7;

      # Tests for the SS-type data
      for (ii, m) in enumerate(N)
        for (jj, n) in enumerate(M)
          for (kk, h) in enumerate(H)

            freq = [10^i for i in range(-6, stop=log10(pi/h), length=1001)]

            AcTrue, BcTrue, CcTrue, DcTrue = rand(m,m), rand(m,n), rand(n,m), rand(n,n)

            ev = eigvals(AcTrue)
            if !isempty(ev[real(ev).<0])
              AcTrue = AcTrue - one(AcTrue) * 2 * maximum(abs.(real(ev[real(ev).<0])))
            end

            Gtrue = ss(AcTrue, BcTrue, CcTrue, DcTrue)
            valA, fA = sigma(Gtrue, freq);
            sysB = hInf_bilinear_s2z(Gtrue, h);
            valB, fB = sigma(sysB, freq)
            sysC = hInf_bilinear_z2s(hInf_bilinear_s2z(Gtrue, h))
            valC, fC = sigma(sysC, freq)
            sysD = hInf_bilinear_s2z(hInf_bilinear_z2s(hInf_bilinear_s2z(Gtrue, h)),h)
            valD, fD = sigma(sysD, freq)

            @test abs(maximum(svd(sysB.B).S)-maximum(svd(sysB.C).S)) < tolerance
            @test abs(maximum(svd(sysC.B).S)-maximum(svd(sysC.C).S)) < tolerance
            @test abs(maximum(svd(sysD.B).S)-maximum(svd(sysD.C).S)) < tolerance

            # Test that the C->D->C and D->C->D results in the same
            @test norm(valA-valC, Inf) < tolerance
            @test norm(valB-valD, Inf) < tolerance
          end
        end
      end
    end
  end
end
