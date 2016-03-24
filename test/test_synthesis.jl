module TestSynthesis
using CustomTest
using ControlSystems

P = tf(1,[1,1])
C = tf([1,1],[1,0])

B = [1]
A = [1,1]
R = [1,1]
S = [1]
T = [1]

@test feedback(P,C) == feedback(P*C) == tf([1,2,1,0],[1,3,3,1,0])
@test feedback2dof(B,A,R,S,T) == tf(B.*T, conv(A,R) + [0;0;conv(B,S)])
@test feedback2dof(P,R,S,T) == tf(B.*T, conv(A,R) + [0;0;conv(B,S)])
end
