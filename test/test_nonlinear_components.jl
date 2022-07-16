using ControlSystems
using ControlSystems: ratelimit

l = 0.5
T = 2
for l in 0.1:0.1:1, T in 1:5
    nl = ratelimit(l; Tf=0.01)
    res = step(nl, T)
    plot(res, plotx=true, plotu=false)
    @test res.y[end] ≈ l*T rtol=0.01
end

##

using ControlSystems
using ControlSystems: deadzone

th = 0.7
@show nl = deadzone(th)
res = lsim(nl, (x,t)->sin(t), 0:0.01:4pi)
# plot(res)
@test minimum(res.y) ≈ -0.3 rtol=1e-3
@test maximum(res.y) ≈ 0.3 rtol=1e-3

@show nl = deadzone(-2th, th)
res = lsim(nl, (x,t)->2sin(t), 0:0.01:4pi)
# plot(res)
@test minimum(res.y) ≈ -0.6 rtol=1e-3
@test maximum(res.y) ≈ 1.3 rtol=1e-3

