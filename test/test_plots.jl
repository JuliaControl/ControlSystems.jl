
module PlotTests
using ControlSystems, Plots
using VisualRegressionTests, ControlExamplePlots
using Test
gr()
default(show=false)

@testset "test_plots" begin
funcs, refs, eps = getexamples()
# Make it easier to pass tests on different systems
# Set to a factor 2 of common errors
eps = [0.15, 0.015, 0.1, 0.01, 0.01, 0.02, 0.01, 0.15, 0.15, 0.01, 0.01]
res = genplots(funcs, refs, eps=eps, popup=false)

##Explicit enumeration for simpler debugging

#"bode.png"
@test  res[1] |> success
#"nyquist.png"
@test  res[2] |> success
#"sigma.png"
@test  res[3] |> success
#"nichols.png"
@test  res[4] |> success
#"step.png"
@test  res[5] |> success
#"impulse.png"
@test  res[6] |> success
#"lsim.png"
@test  res[7] |> success
#"margin.png"
@test  res[8] |> success
#"gangoffour.png"
@test  res[9] |> success
#"pzmap.png"
@test  res[10] |> success
#"rlocus.png"
@test  res[11] |> success
end

end
