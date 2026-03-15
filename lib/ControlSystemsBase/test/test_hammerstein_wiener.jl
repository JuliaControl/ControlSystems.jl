using ControlSystemsBase
using ControlSystemsBase: HammersteinWienerSystem
using SparseArrays
using Test


# @testset "test_nonlinearity_system" begin



# For simplicity, equality of HammersteinWienerSystems are tested over a finite set of frequencies
ω = 0.0:8

@test typeof(promote(nonlinearity(abs2), ss(1))[1]) == HammersteinWienerSystem{Float64}

@test sprint(show, ss(1, 1, 1, 1) * nonlinearity(abs2)) ==
      "HammersteinWienerSystem{Float64}\n\nP: StateSpace{Continuous, Float64}\nA = \n 1.0\nB = \n 0.0  1.0\nC = \n 1.0\n 0.0\nD = \n 0.0  1.0\n 1.0  0.0\n\nContinuous-time state-space model\n\nNonlinearities: Function[abs2]"


P1 = HammersteinWienerSystem(ss(-1.0, 1, 1, 0))
P2 = HammersteinWienerSystem(ss(-2.0, -1, 1, 1)) # (s+1)/(s+2)

# Equality
@test P1 == deepcopy(P1)
@test P1 != deepcopy(P2)


# Random conversions
s = tf("s")
sys1 = HammersteinWienerSystem(1.0 / s)
@test sys1.P.A == sys1.P.D == fill(0, 1, 1)
@test sys1.P.B * sys1.P.C == fill(1, 1, 1)
@test sys1.f == []

## Multiple nonlinearitys
G = ss(1.0) + nonlinearity(abs2) + nonlinearity(abs2)

s11 = feedback(ss(1 / s), nonlinearity(abs2))
s12 = ss(1 / s)
s21 = HammersteinWienerSystem(ss(1 / (s + 1)))
s22 = ss(1 / (s + 10))

s1 = [s11 s12]
s2 = [s21 s22]

f1 = [s1; s2]
f2 = [
    s11 s12
    s21 s22
]
# ↑ these are tested in ControlSystems/test_nonlinear_timeresp


@test propertynames(nonlinearity(abs2)) == (:P, :f, :nu, :ny)



## Test predefined nonlinearities
# saturation
using ControlSystemsBase: saturation
th = 0.7
G = tf(1.0, [1, 1])
@show nl = saturation(th)
Gnl = G * nl
C = tf(1)
Cnl = nl * C
L = G * C
Lnl = G * Cnl


@test (Lnl + 1) == (Lnl + 1.0)
@test (Lnl + 1) isa ControlSystemsBase.HammersteinWienerSystem{Float64}
@test (Lnl - 1) == (Lnl - 1.0)
@test (Lnl - 1) isa ControlSystemsBase.HammersteinWienerSystem{Float64}
@test (1-Lnl) == (1.0 - Lnl)
@test (1-Lnl) isa ControlSystemsBase.HammersteinWienerSystem{Float64}

@test ControlSystemsBase.HammersteinWienerSystem{Float32}(Lnl) isa ControlSystemsBase.HammersteinWienerSystem{Float32}
@test ControlSystemsBase.HammersteinWienerSystem{Float32}(ss(G)) isa ControlSystemsBase.HammersteinWienerSystem{Float32}

n = convert(ControlSystemsBase.HammersteinWienerSystem{Float32}, 1)
@test n isa ControlSystemsBase.HammersteinWienerSystem{Float32}
@test n.D == [1f0;;]
@test n.nx == 0



# offset
using ControlSystemsBase: offset
o = 1.5
G = tf(1.0, [1, 1])
@show nl = offset(o)
# end

w = exp10.(LinRange(-2, 2, 2))
@test all(freqresp(nl, w) .== 1)
@test evalfr(nl, rand())[] == 1

# MIMO offset
o = randn(2)
nl = offset(o)

G = ssrand(2, 1, 1)
nlG = nl * G
@test nlG.ny1 == G.ny

## test linearize
using ControlSystemsBase: linearize
o = 1.5
nl = offset(o)
@test linearize(nl, 1) == ss(1)
@test linearize(nl, 2) == ss(1)

@test linearize(nonlinearity(abs2), 2) == ss(4)
@test linearize(nonlinearity(abs2), 1) == ss(2)
@test linearize(saturation(1), 0) == ss(1)
@test linearize(saturation(1), 3) == ss(0)

@test linearize(saturation(1)tf(1, [1, 1]), 0) == ss(tf(1, [1, 1]))

## Test error on algebraic loop
nl = saturation(0.5)
G = ss(1)


A, B = ControlSystemsBase.linearize((x,u)->x.^2 + sin.(u), [1.0], [2.0])
@test A ≈ [2.0;;]
@test B ≈ [cos(2);;]

## Test nonlinear_components coverage ==========================================
using ControlSystemsBase: Saturation, DeadZone, Offset, Hysteresis, describing_function, deadzone

@testset "nonlinear_components" begin

    # Show methods
    @test sprint(show, Saturation(1.0)) == "saturation(1.0)"
    @test sprint(show, Saturation(-2.0, 1.0)) == "saturation(-2.0, 1.0)"
    @test sprint(show, Offset(1.5)) == "offset(1.5)"
    @test sprint(show, DeadZone(1.0)) == "deadzone(1.0)"
    @test sprint(show, DeadZone(-2.0, 1.0)) == "deadzone(-2.0, 1.0)"
    @test sprint(show, Hysteresis(1.0, 0.5, 20.0)) == "Hysteresis(amplitude=1.0, width=0.5, hardness=20.0)"

    # Hysteresis callable
    h = Hysteresis(1.0, 0.5, 20.0)
    @test isfinite(h(10.0))
    @test h(10.0) ≈ 0.5 * tanh(20.0 * 10.0)
    h_inf = Hysteresis(1.0, 0.5, Inf)
    @test h_inf(10.0) == 0.5
    @test h_inf(-10.0) == -0.5

    # Vector constructors (MIMO)
    @test saturation([1.0, 2.0]) isa HammersteinWienerSystem
    @test deadzone([0.5, 1.0]) isa HammersteinWienerSystem
    @test offset([1.0, 2.0]) isa HammersteinWienerSystem

    # describing_function — numerical
    df_num = describing_function(x -> clamp(x, -1, 1), 2.0)
    df_ana = describing_function(Saturation(1.0), 2.0)
    @test df_num ≈ df_ana rtol=1e-3
    @test_throws ArgumentError describing_function(abs, -1.0)

    # describing_function — Saturation analytical
    @test describing_function(Saturation(2.0), 1.0) == 1.0 + 0im  # A ≤ d
    @test describing_function(Saturation(1.0), 2.0) ≈ 0.6089977810442294
    @test describing_function(Saturation(-0.5, 1.0), 2.0) isa Complex  # asymmetric fallback

    # describing_function — DeadZone analytical
    @test describing_function(DeadZone(2.0), 1.0) == 0.0 + 0im  # A ≤ d
    df_sat = describing_function(Saturation(1.0), 2.0)
    df_dz = describing_function(DeadZone(1.0), 2.0)
    @test df_sat + df_dz ≈ 1.0 + 0im  # complementary
    @test describing_function(DeadZone(-0.5, 1.0), 2.0) isa Complex  # asymmetric fallback

    # describing_function — Hysteresis analytical
    @test describing_function(Hysteresis(1.0, 0.5, 20.0), 0.3) == 0.0 + 0im  # A ≤ width
    h = Hysteresis(1.0, 0.5, 20.0)
    A = 2.0
    r = 0.5 / A
    expected = (4 * 1.0 / (π * A)) * complex(sqrt(1 - r^2), -r)
    @test describing_function(h, A) ≈ expected

    # describing_function — HammersteinWienerSystem
    nl = saturation(1.0)
    @test describing_function(nl, 2.0) ≈ describing_function(Saturation(1.0), 2.0)
    @test_throws ErrorException describing_function(nonlinearity(abs2) + nonlinearity(abs), 1.0)

end
