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
# eps = [0.15, 0.015, 0.1, 0.01, 0.01, 0.02, 0.01, 0.15, 0.15, 0.01, 0.01]
# res = genplots(funcs, refs, eps=eps, popup=false)

# ##Explicit enumeration for simpler debugging
# PROCESSING_ERROR = ControlExamplePlots.PROCESSING_ERROR

# ImageMagick needed to run image comparisons, but it doesnt work on Julia 1.0
# see eg. https://github.com/JuliaIO/ImageMagick.jl/issues/142

#"bode.png"
@test  funcs[1]() isa Plots.Plot
#"nyquist.png"
@test  funcs[2]() isa Plots.Plot
#"sigma.png"
@test  funcs[3]() isa Plots.Plot
#"nichols.png"
@test  funcs[4]() isa Plots.Plot
#"step.png"
@test  funcs[5]() isa Plots.Plot
#"impulse.png"
@test  funcs[6]() isa Plots.Plot
#"lsim.png"
@test  funcs[7]() isa Plots.Plot
#"margin.png"
@test  funcs[8]() isa Plots.Plot
#"gangoffour.png"
@test  funcs[9]() isa Plots.Plot
#"pzmap.png"
@test  funcs[10]() isa Plots.Plot
#"rlocus.png"
@test  funcs[11]() isa Plots.Plot
end

end
