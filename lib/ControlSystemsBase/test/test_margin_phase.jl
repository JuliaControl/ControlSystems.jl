# Number of roots of 1 + G(s)e^{-sθ} with Re s > σ for rational G, computed with the argument principle along the line Re s = σ. This is used as a reference that is independent of rational delay approximations. The frequency grid is refined where the curve changes rapidly relative to its distance from the origin, where the delay rotates the curve by more than half a radian while |G| is not small, and around the imaginary parts of the poles of G.
function n_unstable_delayed(G, θ; σ=1e-9, wmax=1e7)
    g = ss(G)
    pg = poles(g)
    P = count(p -> real(p) > σ, pg)
    f(ω) = 1 + evalfr(g, σ + im*ω)[1]*exp(-(σ + im*ω)*θ)
    w = exp10.(range(-10, log10(wmax), length=2000))
    w = [-reverse(w); 0; w]
    for p in pg
        append!(w, imag(p) .+ (abs(real(p)) + σ) .* range(-50, 50, length=201))
    end
    w = sort!(unique!(filter!(ω -> abs(ω) <= wmax, w)))
    total = 0.0
    stack = [(w[i], w[i+1], f(w[i]), f(w[i+1])) for i in length(w)-1:-1:1]
    while !isempty(stack)
        a, b, fa, fb = pop!(stack)
        m = (a + b)/2
        fm = f(m)
        d = angle(fb/fa)
        coarse = abs(d) > 0.2 || abs(fm - (fa + fb)/2) > 0.1min(abs(fa), abs(fb)) ||
                 (θ*(b - a) > 0.5 && max(abs(fa - 1), abs(fb - 1)) > 0.25)
        if coarse && b - a > 1e-14*max(abs(a), abs(b), 1)
            push!(stack, (m, b, fm, fb))
            push!(stack, (a, m, fa, fm))
        else
            total += d
        end
    end
    round(Int, P - total/2π)
end

# Number of closed-loop poles outside the unit circle for the discrete-time loop G(z)z^{-d}
n_unstable_discrete(G, d) = count(p -> abs(p) > 1 + 1e-9, poles(feedback(ss(G)*ss(tf(1, [1; zeros(d)], G.Ts)))))

# Check every interval returned by delay_stability_intervals against the reference, with the total delay θ = T + τ
function check_delay_intervals(G, T, intervals)
    for (a, b) in intervals
        hi = isfinite(b) ? b : a + 1
        for x in range(a, hi, length=5)[2:end-1]
            @test n_unstable_delayed(G, T + x) == 0
        end
        T + a > 1e-6 && @test n_unstable_delayed(G, T + a - 1e-3*(b - a)) > 0
        isfinite(b) && @test n_unstable_delayed(G, T + b + 1e-3*(b - a)) > 0
    end
    # Delays outside the returned intervals result in an unstable closed loop, so no stable interval is omitted
    bounds = filter(isfinite, [x for iv in intervals for x in iv])
    upper = 1.5*maximum(bounds; init=1.0) + 1
    for x in range(-T, upper, length=41)
        inside = any(a < x < b for (a, b) in intervals)
        near = any(abs(x - y) < 1e-3*(1 + abs(y)) for y in bounds)
        (inside || near) || @test n_unstable_delayed(G, T + x) > 0
    end
end

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
        @test τ ≈ (2atan(3) - π/2 - atan(ε/3))/3 rtol=1e-8
        intervals = delay_stability_intervals(L)
        @test only(intervals) == (0, τ)
    end
    L = 2.7(s+1)^2/((s-1e-4)^2*(s+1e-4))
    check_delay_intervals(L, 0, delay_stability_intervals(L))

    # Starting the grid above the displaced poles selects a different Bode phase branch.
    L = ss(2.7(s+1)^2/((s-1e-4)^2*(s+1e-4)))
    local_grid = [2.0, 3.0, 4.0]
    wide_grid = [0.0; exp10.(range(-7, 0; length=100)); local_grid]
    @test margin(L, local_grid).pm[] ≈ margin(L, wide_grid).pm[]

    # The integrator classification does not change the margins.
    L = (s+0.1)^2/(s^3*(0.01s+1))
    @test margin(L; adjust_phase_start=true).pm[] == margin(L; adjust_phase_start=false).pm[]
    @test margin((s+1e-9)/(s^2*(s+1))).pm[] ≈ 51.8273 atol=1e-2
    @test margin(tf([0.01], [1, -(1-1e-9)], 0.01)).pm[] ≈ 89.7135 atol=1e-2
    @test margin(tf(Float32[2.7, 5.4, 2.7], Float32[1, -1f-4, -1f-8, 1f-12])).pm[] ≈ 53.138 atol=1e-2
    dm32 = delaymargin(tf(Float32[2.7, 5.4, 2.7], Float32[1, -1f-4, -1f-8, 1f-12]))
    @test dm32 isa Float32
    @test dm32 ≈ (2atan(3) - π/2 - atan(1e-4/3))/3 rtol=1e-4

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
    @test delaymargin(Ld) ≈ expected rtol=1e-8
    @test delay_stability_intervals(Ld) == [0:2]
    @test [n_unstable_discrete(Ld, d) == 0 for d in 0:5] == [true, true, true, false, false, false]

    # The smallest phase margin and smallest delay occur at different crossovers.
    L = 0.5/(s^2+0.2s+1) + 20/(s^2+2s+100)
    @test isstable(feedback(L))
    @test length(margin(L; allMargins=true).wpm[]) == 4
    @test margin(L).wpm[] ≈ 1.16742 rtol=1e-3
    τ = delaymargin(L)
    @test n_unstable_delayed(L, 0.995τ) == 0
    @test n_unstable_delayed(L, 1.005τ) > 0
    # The closed loop regains stability for larger delays.
    intervals = delay_stability_intervals(L)
    @test length(intervals) > 1
    @test intervals[1] == (0, τ)
    check_delay_intervals(L, 0, intervals)

    # Select the margin closest to zero, rather than the most negative margin.
    m = margin(-L)
    @test m.pm[] ≈ -9.37304 atol=0.02
    @test m.wpm[] ≈ 0.62356 rtol=1e-3

    @test delaymargin(tf(0.1, [1, 1])) == Inf
    @test delaymargin(tf(1, [1, 1])) == Inf
    @test only(delay_stability_intervals(tf(0.1, [1, 1]))) == (0, Inf)
    # The closed-loop pole at the origin does not depend on the delay.
    @test delaymargin(tf(-1, [1, 1])) == 0
    @test isempty(delay_stability_intervals(tf(-1, [1, 1])))
end

@testset "Delay margin and closed-loop stability" begin
    s = tf("s")

    # Open-loop unstable, and widely separated time scales
    for L in (4(s+1)^2/((s-1)^2*(0.1s+1)), 1e7/((s+1e-4)*(s+1e4)))
        τ = delaymargin(L)
        @test n_unstable_delayed(L, 0.995τ) == 0
        @test n_unstable_delayed(L, 1.005τ) > 0
        @test only(delay_stability_intervals(L)) == (0, τ)
    end
    @test delaymargin(4(s+1)^2/((s-1)^2*(0.1s+1))) ≈ 0.04442 rtol=1e-3
    @test delaymargin(1e7/((s+1e-4)*(s+1e4))) ≈ 0.0014789 rtol=1e-3

    # Unstable closed loops without additional delay do not have a delay margin.
    for L in (27/(s+1)^3, 0.999/(s*(s^2+0.002s+1)), 2(s^2+0.01s+1.0201)/(s*(s^2+0.01s+1)*(s+1)))
        @test n_unstable_delayed(L, 0) > 0
        @test isnan(delaymargin(L))
        @test isempty(delay_stability_intervals(L))
    end

    # An unstable closed loop that is stabilized by an interval of delays
    L = 0.5/(s^2 - 0.1s + 1)
    @test isnan(delaymargin(L))
    intervals = delay_stability_intervals(L)
    @test length(intervals) == 1
    @test intervals[1][1] > 0
    check_delay_intervals(L, 0, intervals)

    # Input or output delays are part of the nominal loop.
    G = 2/(s+1)
    T = 0.3
    for L in (G*delay(T), delay(T)*G, ss(G)*delay(T))
        @test delaymargin(L) ≈ delaymargin(G) - T rtol=1e-8
        @test only(delay_stability_intervals(L)) == (-T, delaymargin(L))
    end
    check_delay_intervals(G, T, delay_stability_intervals(G*delay(T)))

    # The phase margin of an integrator with delay π/2 is zero for unit gain.
    for k in (0.9, 0.999)
        @test delaymargin(k/s*delay(π/2)) ≈ π/(2k) - π/2 rtol=1e-8
    end
    @test delaymargin(1/s*delay(π/2)) == 0
    @test isnan(delaymargin(1.001/s*delay(π/2)))
    @test delay_stability_intervals(1.001/s*delay(π/2))[1][2] < 0

    # A delayed loop of neutral type, |L(j∞)| > 1, is unstable for every positive delay.
    L = (2s+1)/(s+3)
    @test delaymargin(L) == 0
    @test isnan(delaymargin(L*delay(0.1)))
    @test isempty(delay_stability_intervals(L))

    # Internal delays: nominal stability cannot be determined.
    L = feedback(1/(s+1)*delay(0.1))*(5/(s+1))
    @test_logs (:warn, r"cannot be determined") delaymargin(L)
    @test_throws ArgumentError delay_stability_intervals(L)
end

@testset "Delay margin edge cases" begin
    s = tf("s")

    # The reference counts roots of lightly damped and strongly delayed loops.
    @test n_unstable_delayed(tf([4e-5], [1, 2e-5, 1]), 1.0) == 2
    @test n_unstable_delayed(2/(s+1), 1000.0) == 552

    # Nominal closed-loop poles at ±j, with a crossing at zero delay that leaves the right half-plane
    d = 0.5352016707409908
    G = tf([-d, -0.5], [1, d, 1.5])
    @test delaymargin(G) == 0
    intervals = delay_stability_intervals(G)
    @test intervals[1][1] == 0
    @test intervals[1][2] ≈ 3.6172795537 rtol=1e-6
    check_delay_intervals(G, 0, intervals)
    @test n_unstable_delayed(G, 4.0) > 0
    @test isnan(delaymargin(G*delay(4)))

    # Slow unstable closed-loop poles
    @test isnan(delaymargin(1e-9/(s-2e-9)))
    @test isempty(delay_stability_intervals(1e-9/(s-2e-9)))
    @test isnan(delaymargin(tf(Float32[1e-4], Float32[1, -2f-4])))
    @test isnan(delaymargin(tf([1e-10], [1, -(1+3e-10)], 1)))

    # Discrete-time loops with a gain larger than one at the Nyquist frequency
    G = ss(2.0, 1.0)
    @test delaymargin(G) == 1
    @test delay_stability_intervals(G) == [0:0]
    G = tf([-1.2], [1, 0.5], 1)
    @test delaymargin(G) == 1
    @test n_unstable_discrete(G, 1) > 0
    ranges = delay_stability_intervals(G)
    @test [any(d in r for r in ranges) for d in 0:30] == [n_unstable_discrete(G, d) == 0 for d in 0:30]
    # A lightly damped discrete-time loop with stability switches
    G = c2d(ss(0.5/(s^2+0.2s+1)) + ss(20/(s^2+2s+100)), 0.05)
    ranges = delay_stability_intervals(G)
    @test length(ranges) > 1
    @test [any(d in r for r in ranges) for d in 0:150] == [n_unstable_discrete(G, d) == 0 for d in 0:150]
    @test floor(delaymargin(G)) == last(ranges[1])

    # Gain crossover at the Nyquist frequency with a closed-loop pole at z = -1
    @test delaymargin(tf([0.5], [1, 0.5], 1)) == 0
    @test delaymargin(tf([-0.5], [1, 0.5], 1)) == 1

    # Complex coefficients: the crossover is located at a negative frequency
    G = ss(-1-2im, 1, 2, 0)
    τ = delaymargin(G)
    @test τ ≈ (2π/3)/(2+√3) rtol=1e-8
    @test n_unstable_delayed(G, 0.995τ) == 0
    @test n_unstable_delayed(G, 1.005τ) > 0
    check_delay_intervals(G, 0, delay_stability_intervals(G))

    # A closed-loop pole at the origin for every delay, with a computed DC gain that differs from -1 by rounding
    for k in (49.0, 98.0)
        @test delaymargin(tf([-k], [1, k])) == 0
        @test isempty(delay_stability_intervals(tf([-k], [1, k])))
    end

    # Narrow resonances
    @test delaymargin(tf([4e-7], [1, 2e-7, 1])) ≈ π/6 rtol=1e-6
    @test delaymargin(tf(Float32[0.004], Float32[1, 0.002, 1])) ≈ 0.52369187 rtol=1e-3

    # |L(j∞)| = 1: asymptotically, but not exponentially, stable for positive delays
    @test delaymargin((s+1)/(s+2)) == 0
    @test delaymargin(ss(1.0)) == 0
end
