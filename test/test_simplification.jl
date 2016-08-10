module TestSimplification
using CustomTest
using Base.Test
using ControlSystems

## SMINREAL ##
G = ss([-5 0 0 0; 0 -1 -2.5 0; 0 4 0 0; 0 0 0 -6], [2 0; 0 1; 0 0; 0 2],
       [0 3 0 0; -2 0 0 1], [0 0; 1 0])
@test sminreal(G) == G
@test sminreal(G[1, 1]) == ss(0)
@test sminreal(G[1, 2]) == ss([-1 -2.5; 4 0], [1; 0], [3 0], [0])
@test sminreal(G[2, 1]) == ss([-5], [2], [-2], [1])
@test sminreal(G[2, 2]) == ss([-6], [2], [1], [0])

end
