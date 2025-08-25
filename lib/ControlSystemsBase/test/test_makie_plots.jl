using ControlSystemsBase
using Test
using LinearAlgebra
using CairoMakie

@testset "Makie Plot Tests" begin
    # Create test systems
    P = tf([1.2], [1, 2, 1])*tf(1, [10, 1])
    P2 = tf([1, 2], [1, 3, 2])
    Pss = ss(P)
    Pmimo = [P P2; P2 P]
    Pdiscrete = c2d(P, 0.1)
    
    # Test frequency vector
    w = exp10.(range(-2, stop=2, length=100))
    
    @testset "pzmap" begin
        @test_nowarn begin
            fig = CSMakie.pzmap(P)
            CSMakie.pzmap!(fig[1,2], P)
            @test fig isa Makie.Figure
        end
        @test_nowarn begin
            fig = CSMakie.pzmap([P, P2])
            @test fig isa Makie.Figure
        end
        @test_nowarn begin
            fig = CSMakie.pzmap(Pdiscrete; hz=true)
            @test fig isa Makie.Figure
        end
    end
    
    @testset "bodeplot" begin
        @test_nowarn begin
            fig = CSMakie.bodeplot(P)
            CSMakie.bodeplot!(fig[1,2], P)
            @test fig isa Makie.Figure
        end
        @test_nowarn begin
            fig = CSMakie.bodeplot(P, w)
            @test fig isa Makie.Figure
        end
        @test_nowarn begin
            fig = CSMakie.bodeplot([P, P2]; plotphase=false)
            @test fig isa Makie.Figure
        end
        @test_nowarn begin
            fig = CSMakie.bodeplot(Pmimo; hz=true)
            @test fig isa Makie.Figure
        end
        # Test dB scale
        @test_nowarn begin
            ControlSystemsBase.setPlotScale("dB")
            fig = CSMakie.bodeplot(P)
            @test fig isa Makie.Figure
            ControlSystemsBase.setPlotScale("log10")
        end
    end
    
    @testset "nyquistplot" begin
        @test_nowarn begin
            fig = CSMakie.nyquistplot(P)
            CSMakie.nyquistplot!(fig[1,2], P)
            @test fig isa Makie.Figure
        end
        @test_nowarn begin
            fig = CSMakie.nyquistplot(P, w; Ms_circles=[1.5, 2.0])
            @test fig isa Makie.Figure
        end
        @test_nowarn begin
            fig = CSMakie.nyquistplot([P, P2]; unit_circle=true)
            @test fig isa Makie.Figure
        end
    end
    
    @testset "sigmaplot" begin
        @test_nowarn begin
            fig = CSMakie.sigmaplot(P)
            CSMakie.sigmaplot!(fig[1,2], P)
            @test fig isa Makie.Figure
        end
        @test_nowarn begin
            fig = CSMakie.sigmaplot(Pmimo; extrema=true)
            @test fig isa Makie.Figure
        end
    end
    
    @testset "marginplot" begin
        @test_nowarn begin
            fig = CSMakie.marginplot(P)
            CSMakie.marginplot!(fig[1,2], P)
            @test fig isa Makie.Figure
        end
        @test_nowarn begin
            fig = CSMakie.marginplot(P; plotphase=false)
            @test fig isa Makie.Figure
        end
    end
    
    @testset "rlocusplot" begin
        @test_nowarn begin
            fig = CSMakie.rlocusplot(P)
            CSMakie.rlocusplot!(fig[1,2], P)
            @test fig isa Makie.Figure
        end
        @test_nowarn begin
            fig = CSMakie.rlocusplot(P, 100)
            @test fig isa Makie.Figure
        end
    end
    
    @testset "rgaplot" begin
        @test_nowarn begin
            fig = CSMakie.rgaplot(Pmimo)
            CSMakie.rgaplot!(fig[2,1:2], Pmimo)
            @test fig isa Makie.Figure
        end
        @test_nowarn begin
            fig = CSMakie.rgaplot(Pmimo, w; hz=true)
            @test fig isa Makie.Figure
        end
    end
    
    @testset "SimResult plot" begin
        # Create a simple simulation result
        t = 0:0.01:5
        res = step(P, t)
        @test_nowarn begin
            fig = plot(res)
            plot!(fig[1,2], res)
            @test fig isa Makie.Figure
        end
        @test_nowarn begin
            fig = plot(res; plotu=true, plotx=true)
            plot!(fig[1,2], res, plotu=true, plotx=true)
            @test fig isa Makie.Figure
        end
    end
    
    @testset "StepInfo plot" begin
        res = step(P, 100)
        si = stepinfo(res)
        @test_nowarn begin
            fig = plot(si)
            plot!(fig[1,2], si)
            @test fig isa Makie.Figure
        end
    end
    
    @testset "leadlinkcurve" begin
        @test_nowarn begin
            fig = CSMakie.leadlinkcurve()
            CSMakie.leadlinkcurve!(fig[1,2])
            @test fig isa Makie.Figure
        end
        @test_nowarn begin
            fig = CSMakie.leadlinkcurve(2)
            @test fig isa Makie.Figure
        end
    end
    
    @testset "nicholsplot" begin
        @test_nowarn begin
            fig = CSMakie.nicholsplot(P)
            @test fig isa Makie.Figure
        end
    end
end