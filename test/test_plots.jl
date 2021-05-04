# Set plot globals
ENV["PLOTS_TEST"] = "true"
ENV["GKSwstype"] = "nul"

using Plots
gr()
default(show=false)

# This function show mirror that in ControlExamplePlots.jl/genplots.jl
# to make sure that the plots in these tests can be tested for accuracy
"""funcs, names = getexamples()
Get the example functions and names
"""
function getexamples()
    tf1 = tf([1],[1,1])
    tf2 = tf([1/5,2],[1,1,1])
    sys = [tf1 tf2]
    sysss = ss([-1 2; 0 1], [1 0; 1 1], [1 0; 0 1], [0.1 0; 0 -0.2])

    ws = 10.0 .^range(-2,stop=2,length=200)
    ts = 0:0.01:5
    bodegen() = begin
      setPlotScale("dB")
      bodeplot(sys,ws)
    end
    nyquistgen() = nyquistplot(sysss,ws)
    sigmagen() = sigmaplot(sysss,ws)
    #Only siso for now
    nicholsgen() = nicholsplot(tf1,ws)

    stepgen() = stepplot(sys, ts[end], l=(:dash, 4))
    impulsegen() = impulseplot(sys, ts[end], l=:blue)
    L = lqr(sysss.A, sysss.B, [1 0; 0 1], [1 0; 0 1])
    lsimgen() = lsimplot(sysss, (x,i)->-L*x, ts, [1;2])

    margingen() = marginplot([tf1, tf2], ws)
    gangoffourgen() = begin
      setPlotScale("log10");
      gangoffourplot(tf1, [tf(1), tf(5)])
    end
    pzmapgen() = pzmap(tf2, xlims=(-15,5))
    rlocusgen() = rlocusplot(tf2)

    refs = ["bode.png", "nyquist.png", "sigma.png", "nichols.png", "step.png",
            "impulse.png", "lsim.png", "margin.png", "gangoffour.png", "pzmap.png", "rlocus.png"]
    funcs = [bodegen, nyquistgen, sigmagen, nicholsgen, stepgen,
             impulsegen, lsimgen, margingen, gangoffourgen, pzmapgen, rlocusgen]

    funcs, refs
end


@testset "test_plots" begin
  funcs, names = getexamples()

  for i in eachindex(funcs)
    println("run_plot_$(i): $(names[i])")
    @test funcs[i]() isa Plots.Plot
  end

end
