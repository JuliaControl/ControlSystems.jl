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
  end
end
