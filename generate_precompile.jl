G1 = ssrand(1,1,1)
G2 = tf(G1)
G3 = ssrand(1,1,1,Ts=1)
G4 = tf(G3)

u = randn(1,5)
w = exp10.(LinRange(-2, 2, 20))
for G in (G1, G2, G3, G4)
    isdiscrete(G) && lsim(G, u)
    bode(G, w)
    nyquist(G, w)
end
