@testset "test_simplification" begin
## SMINREAL ##
G = ss([-5 0 0 0; 0 -1 -2.5 0; 0 4 0 0; 0 0 0 -6], [2 0; 0 1; 0 0; 0 2],
       [0 3 0 0; -2 0 0 1], [0 0; 1 0])
@test sminreal(G) == G
@test sminreal(G[1, 1]) == ss(0)
@test sminreal(G[1, 2]) == ss([-1 -2.5; 4 0], [1; 0], [3 0], [0])
@test sminreal(G[2, 1]) == ss([-5], [2], [-2], [1])
@test sminreal(G[2, 2]) == ss([-6], [2], [1], [0])


## MINREAL ##

s = tf("s")
P = [1/(s+1) 2/(s+3); 1/(s+1) 1/(s+1)]
sys = ss(P)
sysmin = minreal(sys)

@test size(sysmin.A,1) == 3 # Test that the reduction of sys worked

t = 0:0.1:10
y1,x1 = step(sys,t)[[1,3]]
y2,x2 = step(sysmin,t)[[1,3]]
@test sum(abs2,y1.-y2) < √(eps()) # Test that the output from the two systems are the same

end

@testset "test_promotion" begin
P = tf(1)
s = ss(1)
TP = typeof(P)
Ts = typeof(s)
@test promote_rule(TP,Ts) == Ts

P = tf(1.)
s = ss(1.)
TP = typeof(P)
Ts = typeof(s)
@test promote_rule(Ts,TP) == Ts


P = tf(1.)
s = ss(1)
TP = typeof(P)
Ts = typeof(s)
@test promote_rule(Ts,TP) == typeof(ss(1.))


P = zpk(1)
s = ss(1)
TP = typeof(P)
Ts = typeof(s)
@test promote_rule(Ts,TP) == Ts

P = zpk(1)
P2 = zpk(1.)
TP = typeof(P)
TP2 = typeof(P2)
@test promote_rule(TP,TP2) == TP2

end
