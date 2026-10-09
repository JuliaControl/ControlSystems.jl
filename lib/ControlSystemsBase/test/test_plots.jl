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

@testset "nyquistplot limits" begin
  using ControlSystemsBase: nyquist_limits, nyquist_limit_mask, _nyquist_limit_circles
  s = tf("s")
  w = exp10.(range(-3, 3, length=1000))
  flat(lims) = [l for lims in lims for l in (lims[1]..., lims[2]...)]
  contains(lims, x, y) = lims[1][1] <= x <= lims[1][2] && lims[2][1] <= y <= lims[2][2]

  # The limits contain the entire curve of a system with a high gain
  G = 1000/(s+1)^3
  redata, imdata = nyquist(G, w)
  lims = nyquist_limits(G, w)[1, 1]
  @test all(contains.(Ref(lims), redata[1, 1, :], imdata[1, 1, :]))
  @test contains(lims, -1.5, 0.5) && contains(lims, -0.5, -0.5)

  # A curve that is located far from the critical point
  G = ss(100) + tf(1, [1, 1])
  lims = nyquist_limits(G, w)[1, 1]
  @test contains(lims, 101 - 1e-3, 0) && contains(lims, -1.5, 0)
  @test nyquistplot(G) isa Plots.Plot

  # The limits of each channel are computed separately, and contain the curves of all systems
  Gmimo = [tf(1, [1, 1]) tf(10, [1, 1]); tf(0.1, [1, 1]) tf(10, [1, 2, 1])]
  lims = nyquist_limits(Gmimo, w)
  @test lims[1, 1][1][2] < 2
  @test lims[1, 2][1][2] > 10
  @test lims[2, 1][1][2] < 0.5 && lims[2, 1] != lims[1, 1]
  @test nyquist_limits([tf(1, [1, 1]), tf(10, [1, 1])], w)[1, 1] == lims[1, 2]

  # Invariance under a change of the time unit
  for G in (ss(1/(s*(s+1))), ss(1/((s^2 + 1)*(s+1))), ss(1000/(s+1)^3))
    @test flat(nyquist_limits(G, w)) ≈ flat(nyquist_limits(ControlSystemsBase.time_scale(G, 100), 100 .* w))
  end

  # The frequencies at which the factor contributed by the poles on the imaginary axis exceeds max_factor are excluded
  @test nyquist_limit_mask(1/(s*(s+1)), w) == (w .>= 1/5)
  @test nyquist_limit_mask(1/(s^2*(s+1)), w) == (w .>= 1/sqrt(5))
  @test nyquist_limit_mask((s+1)/(s*(s+10)), w; max_factor=4) == (w .>= 1/4)
  @test nyquist_limit_mask(1/((s^2 + 1)*(s+1)), w) == (abs.(1 .- w.^2) .>= 1/5)
  @test nyquist_limit_mask(c2d(ss(1/(s*(s+1))), 0.01), w[w .< 300]) == (w[w .< 300] .>= 1/5)
  @test nyquist_limit_mask(1/(s+1)^3, w) == trues(length(w))
  @test !any(nyquist_limit_mask(1/s^2, w)) # No intrinsic scale

  # A system with an integrator has finite limits that contain the critical point
  lims = nyquist_limits(1/(s*(s+1)), w)[1, 1]
  @test all(isfinite, (lims[1]..., lims[2]...))
  @test contains(lims, -1, 0)
  @test -6 < lims[2][1] < imag(1/(0.2im*(0.2im+1)))
  lims = nyquist_limits(1/s^2, w)[1, 1]
  @test flat([lims]) ≈ [-1.65, 0.15, -0.6, 0.6]
  @test nyquist_limits(delay(1)*tf(1, [1, 1, 0]), w)[1, 1][2][1] > -10

  @test nyquist_limits(tf(0.1, [1, 1]), w; critical_point=-2)[1, 1][1][1] < -3
  @test length(_nyquist_limit_circles([1.2], [1.05, 2], Float64[], true)) == 3 # The Mt circle of radius > 2 is excluded
  lims = nyquist_limits(tf(0.1, [1, 1]), w; circles=_nyquist_limit_circles([], [2], [], false))[1, 1]
  @test contains(lims, -2, 0)

  # The Plots recipe sets the limits of each subplot, user-provided limits take precedence
  p = nyquistplot(Gmimo, w)
  lims = nyquist_limits(Gmimo, w)
  for i = 1:2, j = 1:2
    sp = p.subplots[LinearIndices((2, 2))[j, i]]
    @test Plots.xlims(sp) == lims[i, j][1]
    @test Plots.ylims(sp) == lims[i, j][2]
  end
  p = nyquistplot(1/(s*(s+1)), xlims=(-3, 3))
  @test Plots.xlims(p.subplots[1]) == (-3, 3)

  # The default limits are not computed for number types for which the poles are not available (the poles of a BigFloat system require GenericSchur)
  Gbig = tf(big(1.0), big.([1.0, 2, 1]))
  @test !ControlSystemsBase._nyquist_limits_available([Gbig])
  @test ControlSystemsBase._nyquist_limits_available([tf(1, [1, 1]), ss(1.0f0)])
  @test nyquistplot(Gbig, w) isa Plots.Plot
  p = nyquistplot(Gbig, w, xlims=(-3, 3), ylims=(-2, 2))
  @test Plots.xlims(p.subplots[1]) == (-3, 3) && Plots.ylims(p.subplots[1]) == (-2, 2)
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
