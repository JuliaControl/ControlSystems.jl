@testset "Phase margin branches" begin
    s = tf("s")
    for ε in (0.0, 1e-6, 1e-4), representation in (tf, zpk, ss)
        L = representation(2.7(s+1)^2/((s-ε)^2*(s+ε)))
        @test isstable(feedback(L))
        ε == 0 || @test ControlSystemsBase.count_integrators(L) == 0
        for adjust in (false, true)
            m = margin(L; full=true, adjust_phase_start=adjust)
            ω = m.wpm[]
            expected = rad2deg(2atan(ω) - π/2 - atan(ε/ω))
            @test m.pm[] ≈ expected atol=1e-7
            @test margin(L; adjust_phase_start=adjust).pm[] ≈ m.pm[]
            all = margin(L; allMargins=true, adjust_phase_start=adjust)
            @test only(all.pm[]) ≈ m.pm[]
            @test only(all.wpm[]) ≈ ω
            @test rem(m.fullPhase[] - m.pm[] + 180, 360, RoundNearest) ≈ 0 atol=0.02
        end
        τ = delaymargin(L)
        @test τ ≈ (2atan(3) - π/2 - atan(ε/3))/3 rtol=5e-4
        @test isstable(feedback(L*pade(0.95τ, 6)))
        @test !isstable(feedback(L*pade(1.05τ, 6)))
    end

    # Starting the grid above the displaced poles selects a different Bode phase branch.
    L = ss(2.7(s+1)^2/((s-1e-4)^2*(s+1e-4)))
    local_grid = [2.0, 3.0, 4.0]
    wide_grid = [0.0; exp10.(range(-7, 0; length=100)); local_grid]
    @test margin(L, local_grid).pm[] ≈ margin(L, wide_grid).pm[]

    # Excess phase lag must produce a negative signed margin.
    L = 27/(s+1)^3
    m = margin(L, [1.0, sqrt(8), 4.0])
    @test m.pm[] ≈ 180 - 3rad2deg(atan(sqrt(8))) atol=0.02
    @test m.pm[] < 0

    # Delay contributes phase revolutions without changing gain crossovers.
    L = 2/(s+1)
    ω = margin(L).wpm[]
    m = margin(L*delay(2π/ω))
    @test m.pm[] ≈ margin(L).pm[] atol=0.02

    # Discrete-time delay margins are measured in sample intervals.
    Ts = 0.2
    Ld = tf([0.5], [1, -1], Ts)
    ω = 2asin(0.25)/Ts
    expected = (π/2 - ω*Ts/2)/(ω*Ts)
    @test delaymargin(Ld) ≈ expected rtol=1e-3

    # The smallest phase margin and smallest delay occur at different crossovers.
    L = 0.5/(s^2+0.2s+1) + 20/(s^2+2s+100)
    @test isstable(feedback(L))
    @test length(margin(L; allMargins=true).wpm[]) == 4
    @test margin(L).wpm[] ≈ 1.16742 rtol=1e-3
    τ = delaymargin(L)
    @test τ ≈ 0.15874 rtol=1e-3
    @test isstable(feedback(L*pade(0.95τ, 6)))
    @test !isstable(feedback(L*pade(1.05τ, 6)))

    # Select the margin closest to zero, rather than the most negative margin.
    m = margin(-L)
    @test m.pm[] ≈ -9.37304 atol=0.02
    @test m.wpm[] ≈ 0.62356 rtol=1e-3

    @test delaymargin(tf(0.1, [1, 1])) == Inf
    @test delaymargin(tf(1, [1, 1])) == Inf
    @test delaymargin(tf(-1, [1, 1])) == 0
end
