module TestLqg
using CustomTest
using ControlSystems


w = logspace(-5,5,1000)
s = tf("s")
P = [1/(s+1) 2/(s+3); 1/(s+1) 1/(s+1)]
sys = ss(P)
sys.B

q = 0
Q1 = 100eye(4)
Q2 = 1eye(2)
R1 = 100eye(4)
R2 = 1eye(2)
Ginit = LQG(sys, Q1, Q2, R1, R2, qQ=q, qR=q)
G = lqg(Ginit)


qQ = 1
qR = 1
Q1 = 1000eye(4)
Q2 = 1eye(2)
R1 = 1eye(6)
R2 = 1eye(2)
Ginit = LQG(sys, Q1, Q2, R1, R2, qQ=qQ, qR=qR, integrator=true)
Gi = lqg(Ginit)

Gcl = G[:cl]
T = G[:T]
S = G[:S]
CS = G[:CS]
PS = G[:PS]


Gcl = Gi[:cl]
T = Gi[:T]
S = Gi[:S]
sigmaplot(S,w)
Plots.gui()

Plots.plotlyjs()
gangoffourplot(Gi)
Plots.gui()


Plots.plot(randn(10))
Plots.gui()


P = 1/(s^2+0.2s+1)
sys = ss(P)

f = begin
    qQ = 1000
    qR = 100
    Q1 = 1eye(2)
    Q2 = 1eye(1)
    R1 = 10000eye(2)
    R2 = 1eye(1)
    Ginit = LQG(sys, Q1, Q2, R1, R2, qQ=qQ, qR=qR, integrator=false)
    G = lqg(Ginit)

    Gcl = G[:cl]
    T = G[:T]
    S = G[:S]
    sigmaplot([S, T],w)
    Plots.gui()
end





end
