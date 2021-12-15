# Example from the paper
using ControlSystems

L = tf(25, [1,10,10,10])
dm = diskmargin(L, 0)
@test dm.ω0 ≈ 1.94   atol=0.02
@test dm.γmin ≈ 0.63    atol=0.02
@test dm.γmax ≈ 1.59    atol=0.02
@test dm.α ≈ 0.46   atol=0.02
plot(dm)
nyquistplot(L)
plot!(dm, nyquist=true)
plot!(Disk(dm))




## Frequency-dependent margin
w = exp10.(LinRange(-2, 2, 500))
dms = diskmargin(L, 0, w)
plot(w, dms)

##
s = tf("s")
L = 6.25*(s + 3)*(s + 5) / (s*(s + 1)^2 *(s^2 + 0.18s + 100))

## αmax > 2
dm = diskmargin(L, 0, 100)
@test dm.γmax == Inf
@test dm.γmin == 0
@test dm.ϕm == 90


w = exp10.(LinRange(-1, 2, 500))
dms = diskmargin(L, 0, w)
plot(w, dms) 



# using IntervalArithmetic
# δ(a=1) = Complex(-a..a, -a..a)
# Δ(n, a) = diagm([δ(a) for _ in 1:n])
# M = [0 1; -0.1 -0.1]
# D = Δ(2, 1)
# 0 ∈ det(I-M*D)