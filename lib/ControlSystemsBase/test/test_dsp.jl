@testset "DSP interoperability" begin
    @info "Testing DSP interoperability"
    import DSP 
    G = DemoSystems.resonant()*DemoSystems.resonant(ω0=2) |> tf
    Gd0 = c2d(DemoSystems.resonant(), 0.1) |> tf
    Gd = c2d(G, 0.1)
    Gs, k = ControlSystemsBase.seriesform(Gd)
    @test k*prod(Gs) ≈ Gd atol=1e-3
    Gds = DSP.SecondOrderSections(Gd)

    u = randn(100)
    uf = DSP.filt(DSP.DF2TFilter(Gds, zeros(2,2)), u)
    uls = lsim(Gd, u').y'
    @test uf[1:end-1] ≈ uls[2:end]

    z,p,k = randn(3), randn(3), randn()

    f = DSP.ZeroPoleGain(z,p,k)
    fcs = zpk(f, 1)

    @test fcs.matrix[1].z ≈ z
    @test fcs.matrix[1].p ≈ p
    @test fcs.matrix[1].k ≈ k

    u = randn(10)
    uf = DSP.filt(f, u)
    uls = lsim(fcs, u').y'
    @test uf ≈ uls






    u = randn(100)
    res = lsim(Gd, u')
    yf = res.y
    @test yf' ≈ DSP.filt(Gd, u)

    # The filtfilt function has a complicated initialization of the filter state, making it hard to compare to twice lsim
    @test length(DSP.filtfilt(Gd, u)) == length(u)
    
end
