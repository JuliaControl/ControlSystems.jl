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

