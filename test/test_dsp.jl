@testset "DSP interoperability" begin
    @info "Testing DSP interoperability"
    import DSP 
    G = DemoSystems.resonant()*DemoSystems.resonant(ω0=2) |> tf
    Gd = c2d(G, 0.1)
    Gs, k = seriesform(Gd)
    @test k*prod(Gs) ≈ Gd
    Gds = DSP.SecondOrderSections(Gd)
    
end