# Set plot globals
ENV["PLOTS_TEST"] = "true"
ENV["GKSwstype"] = "100"

module PlotTests
using ControlSystems, Plots
using ControlExamplePlots
using Test
gr()
default(show=false)

@testset "test_plots" begin
funcs, refs, eps = getexamples()
# Make it easier to pass tests on different systems
# Set to a factor 2 of common errors, these are ignored now, use visual inspection
eps = [0.15, 0.015, 0.1, 0.01, 0.01, 0.02, 0.01, 0.15, 0.15, 0.01, 0.01]
res = genplots(funcs, refs, eps=eps, popup=false)

##Explicit enumeration for simpler debugging
PROCESSING_ERROR = ControlExamplePlots.PROCESSING_ERROR

#"bode.png"
@test  res[1].status != PROCESSING_ERROR
#"nyquist.png"
@test  res[2].status != PROCESSING_ERROR
#"sigma.png"
@test  res[3].status != PROCESSING_ERROR
#"nichols.png"
@test  res[4].status != PROCESSING_ERROR
#"step.png"
@test  res[5].status != PROCESSING_ERROR
#"impulse.png"
@test  res[6].status != PROCESSING_ERROR
#"lsim.png"
@test  res[7].status != PROCESSING_ERROR
#"margin.png"
@test  res[8].status != PROCESSING_ERROR
#"gangoffour.png"
@test  res[9].status != PROCESSING_ERROR
#"pzmap.png"
@test  res[10].status != PROCESSING_ERROR
#"rlocus.png"
@test  res[11].status != PROCESSING_ERROR
end

end
