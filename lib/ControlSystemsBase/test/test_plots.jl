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

    sysd = c2d(ss(sys), 0.01)
    sysssd = c2d(sysss, 0.01)

    ws = 10.0 .^range(-2,stop=2,length=200)
    ts = 0:0.01:5
    bodegen() = begin
      setPlotScale("dB")
      bodeplot(sys,ws)
    end
    nyquistgen() = nyquistplot(sysss,ws, Ms_circles=1.2, Mt_circles=1.2)
    sigmagen() = sigmaplot(sysss,ws)
    #Only siso for now
    nicholsgen() = nicholsplot(tf1,ws)

    stepgen() = plot(step(sysd, ts[end]), l=(:dash, 4))
    impulsegen() = plot(impulse(sysd, ts[end]), l=:blue)
    L = lqr(sysss, [1 0; 0 1], [1 0; 0 1])
    lsimgen() = plot(lsim(sysssd, (x,i)->-L*x, ts; x0=[1;2]), plotu=true)
    plot(lsim.([sysssd, sysssd], (x,i)->-L*x, Ref(ts); x0=[1;2]), plotu=true, plotx=true)

    margingen() = marginplot([tf1, tf2], ws)
    gangoffourgen() = begin
      setPlotScale("log10");
      gangoffourplot(tf1, [tf(1), tf(5)])
    end
    pzmapgen() = pzmap(c2d(tf2, 0.1))

    refs = ["bode.png", "nyquist.png", "sigma.png", "nichols.png", "step.png",
            "impulse.png", "lsim.png", "margin.png", "gangoffour.png", "pzmap.png"]
    funcs = [bodegen, nyquistgen, sigmagen, nicholsgen, stepgen,
             impulsegen, lsimgen, margingen, gangoffourgen, pzmapgen]

    funcs, refs
end


@testset "test_plots" begin
  funcs, names = getexamples()

  for i in eachindex(funcs)
    println("run_plot_$(i): $(names[i])")
    @test funcs[i]() isa Plots.Plot
  end

  sys = ssrand(3,3,3)
  sigmaplot(sys, extrema=true)

  Gmimo = ssrand(2,2,2,Ts=1)
  @test_nowarn plot(step(Gmimo, 10), plotx=true)
  # plot!(step(Gmimo[:, 1], 10), plotx=true) # Verify that this plots the same as the "from u(1)" series above

end
