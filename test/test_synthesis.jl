module TestSynthesis
using CustomTest
using ControlSystems

P = tf(1,[1,1])
C = tf([1,1],[1,0])
L = P*C
Lsys = ss(L)


B = [1]
A = [1,1]
R = [1,1]
S = [1]
T = [1]

@test minreal(feedback(P,C),1e-5) == minreal(feedback(L),1e-5) == minreal(tf([1,2,1,0],[1,3,3,1,0]),1e-5)
@test feedback2dof(B,A,R,S,T) == tf(B.*T, conv(A,R) + [0;0;conv(B,S)])
@test feedback2dof(P,R,S,T) == tf(B.*T, conv(A,R) + [0;0;conv(B,S)])
@test isapprox(pole(minreal(tf(feedback(Lsys)),1e-5)) , pole(minreal(feedback(L),1e-5)), atol=1e-5) # This test is ugly, but numerical stability is poor for minreal

end
