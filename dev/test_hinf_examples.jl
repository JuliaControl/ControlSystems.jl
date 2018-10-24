using ControlSystems
using Test
using LinearAlgebra

@testset "test_hinf_examples" begin

  ##############################################################################
  # Test the DC motor example
  ##############################################################################

  # Make sure that the code runs
  @test isa(include("hinf_example_DC.jl"), Nothing)
  f = [10^i for i in range(-7, stop=7, length=201)]
  # Check that the optimal gain is correct
  @test abs(gamma - 4.466640) < 1e-2
  # Check that the closed loop tf satisfies ||Pcl||_∞ < gamma
  valPcl  = sigma(Pcl, f)[1];
  @test all(valPcl .< (gamma+1e-4))
  # Check that ||S/WS||_∞ < gamma
  if isa(WS, LTISystem) || isa(WS, Number)
    valSWS = sigma(S * WS , f)[1];
    @test all(valSWS .< (gamma+1e-4))
  end
  # Check that ||S/WS||_∞ < gamma
  if isa(WU, LTISystem) || isa(WU, Number)
    valKSWU = sigma(KS * WU , f)[1];
    @test all(valKSWU .< (gamma+1e-4))
  end
  # Check that ||S/WS||_∞ < gamma
  if isa(WT, LTISystem) || isa(WT, Number)
    valTWT = sigma(T * WT , f)[1];
    @test all(valTWT .< (gamma+1e-4))
  end

  ################################################################################
  # Test the MIT open courseware example
  ################################################################################

  # Make sure that the code runs
  @test isa(include("hinf_example_MIT.jl"), Nothing)
  f = [10^i for i in range(-7, stop=7, length=201)]
  # Check that the optimal gain is correct
  @test abs(gamma - 1.636084) < 1e-2
  # Check that the closed loop tf satisfies ||Pcl||_∞ < gamma
  valPcl  = sigma(Pcl, f)[1];
  @test all(valPcl .< (gamma+1e-4))
  # Check that ||S/WS||_∞ < gamma
  if isa(WS, LTISystem) || isa(WS, Number)
    valSWS = sigma(S * WS , f)[1];
    @test all(valSWS .< (gamma+1e-4))
  end
  # Check that ||S/WS||_∞ < gamma
  if isa(WU, LTISystem) || isa(WU, Number)
    valKSWU = sigma(KS * WU , f)[1];
    @test all(valKSWU .< (gamma+1e-4))
  end
  # Check that ||S/WS||_∞ < gamma
  if isa(WT, LTISystem) || isa(WT, Number)
    valTWT = sigma(T * WT , f)[1];
    @test all(valTWT .< (gamma+1e-4))
  end


  ################################################################################
  # Test the quad tank example
  ################################################################################

  # Make sure that the code runs
  @test isa(include("hinf_example_tank.jl"), Nothing)
  f = [10^i for i in range(-7, stop=7, length=201)]
  # Check that the optimal gain is correct
  @test abs(gamma - 5.052886) < 1e-2
  # Check that the closed loop tf satisfies ||Pcl||_∞ < gamma
  valPcl  = sigma(Pcl, f)[1];
  @test all(valPcl[:,1] .< (gamma+1e-4))
  # Check that ||S/WS||_∞ < gamma
  if isa(WS, LTISystem) || isa(WS, Number)
    valSWS = sigma(S * WS , f)[1];
    @test all(valSWS[:,1] .< (gamma+1e-4))
  end
  # Check that ||S/WS||_∞ < gamma
  if isa(WU, LTISystem) || isa(WU, Number)
    valKSWU = sigma(KS * WU , f)[1];
    @test all(valKSWU[:,1] .< (gamma+1e-4))
  end
  # Check that ||S/WS||_∞ < gamma
  if isa(WT, LTISystem) || isa(WT, Number)
    valTWT = sigma(T * WT , f)[1];
    @test all(valTWT[:,1] .< (gamma+1e-4))
  end

end
