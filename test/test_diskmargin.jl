# Example from the paper
using ControlSystems
using ControlSystems: diskmargin, γϕcurve

L = tf(25, [1,10,10,10])
dm = diskmargin(L, 0)
@test dm.ω0 ≈ 1.94   atol=0.02
@test dm.γmin ≈ 0.63    atol=0.02
@test dm.γmax ≈ 1.59    atol=0.02
@test dm.α ≈ 0.46   atol=0.02
plot(dm)
