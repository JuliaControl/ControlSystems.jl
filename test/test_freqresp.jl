module TestFreqResp
using CustomTest
using Base.Test
using ControlSystems

## EVALFR ##
H = [tf(0) tf([3, 0],[1, 1, 10]) ; tf([1, 1],[1, 5]) tf([2],[1, 6])]
G = ss([-5 0 0 0; 0 -1 -2.5 0; 0 4 0 0; 0 0 0 -6], [2 0; 0 1; 0 0; 0 2],
       [0 3 0 0; -2 0 0 1], [0 0; 1 0])

@test evalfr(H, -6) == [0.0 -0.45; 5.0 Inf]
@test evalfr(H, -5) == [0.0 -0.5; Inf 2.0]
@test evalfr(H, -1) == [0.0 -0.3; 0.0 0.4]
@test_approx_eq evalfr(H, 0) [0.0 0.0; 0.2 1/3]

@test evalfr(G, -6) == [Inf Inf; Inf Inf]
@test evalfr(G, -5) == [Inf Inf; Inf Inf]
@test_approx_eq evalfr(G, -1) [0.0 -0.3; 0.0 0.4]
@test_approx_eq evalfr(G, 0) [0.0 0.0; 0.2 1/3]

## Shortcut notation for evalfr ##
F = tf([1],[1,0.5],-1)
omega = 2
z = 0.5(1+im)
@test F(omega,true)[1] == 1/(exp(-im*2)+0.5)
@test F(z,false)[1] == 1/(z+0.5)
@test_throws ErrorException F(z,true)

## Test bode, nyquist and sigma
sys = [tf([1,-1], [1,1,1]) 0; 0 tf([1],[1,1])]
f(s) = [(s-1)./(s.^2+s+1) 0; 0 1./(1+s)]
ws = logspace(-2,2,50)
resp = Array(Complex128,50,2,2)
for (i,w) in enumerate(ws)
    resp[i,:,:] = f(im*w)
end

@test bode(sys, ws)[1:2] == (abs(resp), rad2deg(angle(resp)))
@test nyquist(sys, ws)[1:2] == (real(resp), imag(resp))
sigs = Array(Float64,50,2)
for i in eachindex(ws)
    sigs[i,:] =  svdvals(resp[i,:,:])
end
@test sigma(sys, ws)[1] == sigs

#Test default freq vector contains at least half a decade more than all features
p, z = 100, 1/1000
sys2 = [tf([1],[1/p,1]) tf([1/z, 2/z, 1],[1])]
mag, phase, ws2 = bode(sys2)
@test maximum(ws2) >= 2max(p,z)
@test minimum(ws2) <= 0.5min(p,z)

#Test default freq vector not too small
p, z = 1, 1.2
sys2 = [tf([1],[1/p,1]) tf([1/z, 2/z, 1],[1])]
mag, mag, ws2 = bode(sys2)
@test maximum(ws2) >= 5max(p,z)
@test minimum(ws2) <= 0.2min(p,z)
@test length(ws2) > 100
end
