using ControlSystemsBase
using ControlSystemsBase: HammersteinWienerSystem
using SparseArrays


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


