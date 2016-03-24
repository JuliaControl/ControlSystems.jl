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

@test isapprox(minreal(feedback(P,C),1e-5), minreal(feedback(L),1e-5))
@test isapprox(numpoly(minreal(feedback(L),1e-5))[1].a, numpoly(tf(1,[1,1]))[1].a)# This test is ugly, but numerical stability is poor for minreal
@test feedback2dof(B,A,R,S,T) == tf(B.*T, conv(A,R) + [0;0;conv(B,S)])
@test feedback2dof(P,R,S,T) == tf(B.*T, conv(A,R) + [0;0;conv(B,S)])
@test isapprox(pole(minreal(tf(feedback(Lsys)),1e-5)) , pole(minreal(feedback(L),1e-5)), atol=1e-5) 

@test_err feedback(ss(1),ss(1))
@test_err feedback(ss(eye(2), ones(2,2), ones(1,2),0))

end
