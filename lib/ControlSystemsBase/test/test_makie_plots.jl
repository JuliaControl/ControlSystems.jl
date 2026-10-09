using ControlSystemsBase
using Test
using LinearAlgebra
import CairoMakie
import CairoMakie.Makie

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
            fig = CSMakie.nyquistplot(P, w; Ms_circles=[1.5, 2.0], Mt_circles=[1.5, 2.0])
            @test fig isa Makie.Figure
        end
        @test_nowarn begin
            fig = CSMakie.nyquistplot([P, P2]; unit_circle=true)
            @test fig isa Makie.Figure
        end
        @test_nowarn begin
            fig = CSMakie.nyquistplot(P; polar=true, rlimits=(0, 3),
                                      Ms_circles=[1.5], Mt_circles=[1.5], unit_circle=true)
            @test fig isa Makie.Figure
            @test any(x -> x isa Makie.PolarAxis, fig.content)
        end
        @test_nowarn begin
            fig = CSMakie.nyquistplot(Pmimo; polar=true)
            @test fig isa Makie.Figure
        end
        # Extra keyword arguments are forwarded to the axis constructor
        @test_nowarn begin
            fig = CSMakie.nyquistplot(P; polar=true, rlimits=(0, 3),
                                      rticks=0:0.5:3,
                                      thetaminorticks=Makie.IntervalsBetween(3),
                                      thetaminorgridvisible=true, rminorgridvisible=true)
            @test fig isa Makie.Figure
        end
        @test_nowarn begin
            fig = CSMakie.nyquistplot(P; xgridvisible=false, title="Custom title")
            @test fig isa Makie.Figure
        end
        # Unsupported keyword arguments now error instead of being silently ignored
        @test_throws Exception CSMakie.nyquistplot(P; not_a_keyword=true)
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
        # Issue #771 and the numerical cases from the review of #1076.
        G = tf([1.0], [1.0, 13, 40, 0])
        Tss = feedback(ss(tf(1, [1, 1, 0])) * ss(pid(1.0, 1.0, 0.1; Tf=0.01)))
        L = ss(tf([1, 1], [1, 10])) * ss(10.0)
        for sys in (G, Tss, L, DemoSystems.double_mass_model())
            @test_nowarn begin
                fig = CSMakie.marginplot(sys)
                @test fig isa Makie.Figure
            end
        end
        # Explicit ranges and fallback phase branches use the same helper
        # as Plots; verify that both sampling and display options work.
        for sys in (G, tf(0.1, [1.0, 1]), tf(-0.1, [1.0, 1]), tf(1e-6, [1.0, 0, 0, 0, 0, 0]))
            for adaptive in (false, true), hz in (false, true)
                @test_nowarn begin
                    fig = CSMakie.marginplot(sys, w; adaptive, hz)
                    @test fig isa Makie.Figure
                end
            end
        end
        @test_nowarn CSMakie.marginplot(Tss, exp10.(range(-8, 4; length=1000)))
        # Each channel of each system uses its own default frequency vector, see `margin`
        let s = tf("s")
            sys = append(ss(1/((s + 1e-6)*(s + 1))), ss(1e4/(s + 1e4)))
            for systems in (sys, [sys, 2sys])
                @test_nowarn begin
                    fig = CSMakie.marginplot(systems)
                    @test fig isa Makie.Figure
                end
            end
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
            fig = Makie.plot(res)
            Makie.plot!(fig[1,2], res)
            @test fig isa Makie.Figure
        end
        @test_nowarn begin
            fig = Makie.plot(res; plotu=true, plotx=true)
            Makie.plot!(fig[1,2], res, plotu=true, plotx=true)
            @test fig isa Makie.Figure
        end

        res = step(Pmimo, t)
        @test_nowarn begin
            fig = Makie.plot(res; plotu=true)
            @test fig isa Makie.Figure
        end
    end
    
    @testset "StepInfo plot" begin
        res = step(P, 100)
        si = stepinfo(res)
        @test_nowarn begin
            fig = Makie.plot(si)
            Makie.plot!(fig[1,2], si)
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