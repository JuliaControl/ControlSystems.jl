module TestFreqResp
using CustomTest
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
end
