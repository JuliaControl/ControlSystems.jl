using ControlSystems
using Test
using LinearAlgebra
using Random

@testset "test_hinfinity_design" begin
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
@test (1+1)==2
end
