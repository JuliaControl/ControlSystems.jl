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
    nyquistgen() = nyquistplot(sysss,ws, Ms_circles=1.2, Mt_circles=1.2, disk_margin_circles=1.2)
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
  bodeplot(sys, adaptive=false)
  nyquistplot(sys, adaptive=false)

  Gmimo = ssrand(2,2,2,Ts=1)
  @test_nowarn plot(step(Gmimo, 10), plotx=true)
  # plot!(step(Gmimo[:, 1], 10), plotx=true) # Verify that this plots the same as the "from u(1)" series above

  setPlotScale("log10")
  funcs[8]() # test marginlpot with log10 scale
  
  # Test marginplot with hz=true
  @test_nowarn marginplot(tf([2], [1,1])^3, 10.0 .^range(-2,stop=2,length=50); hz=true)
end

@testset "marginplot regressions" begin
  G = tf([1.0], [1.0, 13, 40, 0])
  @test_nowarn marginplot(G; xticks=exp10.(-2:0.5:4))
  Tss = feedback(ss(tf(1, [1, 1, 0])) * ss(pid(1.0, 1.0, 0.1; Tf=0.01)))
  L = ss(tf([1, 1], [1, 10])) * ss(10.0)
  for sys in (Tss, L, DemoSystems.double_mass_model())
    @test_nowarn marginplot(sys)
  end

  # Smoke tests for absent phase margins on different phase branches.
  w = exp10.(range(-1, 3; length=200))
  for sys in (tf(0.1, [1.0, 1]), tf(-0.1, [1.0, 1]), tf(1e-6, [1.0, 0, 0, 0, 0, 0]))
    for adaptive in (false, true), hz in (false, true)
      @test_nowarn marginplot(sys, w; adaptive, hz)
    end
  end

  # Keep all phase-margin arrays aligned when the display is truncated.
  peaks = sum(tf(0.1ω^2, [1, 0.01ω, ω^2]) for ω in (1, 10, 100))
  dense = exp10.(range(-2, 3; length=2000))
  @test length(margin(peaks, dense; allMargins=true).pm[]) == 6
  @test_logs (:warn, r"Only showing .* 5 out of 6 phase margins") marginplot(peaks, dense)
  @test marginplot(Tss, exp10.(range(-8, 4; length=1000))) isa Plots.Plot

  # Each channel of each system uses its own default frequency vector, so the slow pole of channel (1, 1) does not produce crossovers caused by rounding error in channel (2, 2), whose gain equals one at low frequencies
  let s = tf("s")
    sys = append(ss(1/((s + 1e-6)*(s + 1))), ss(1e4/(s + 1e4)))
    @test_nowarn marginplot(sys)
    @test_nowarn marginplot([sys, 2sys])
    setPlotScale("dB")
    try
      @test_nowarn marginplot(sys)
    finally
      setPlotScale("log10")
    end
  end
end
