using Plots

# Frequency domain analysis of a simple system with time delays

s = tf("s")

P = delay(0.2) * ss(1/(s+1))

K = 3; Ti = 0.3;
C = DelayLtiSystem(ss(K*(1 + 1/s)))

ω = exp10.(LinRange(-2,2,500))

L_fr = freqresp(C*P, ω)[:]
plot(real(L_fr), imag(L_fr), xlim=[-2,1], ylim=[-2,2], title="Nyquist curve")

G_yd = feedback(P, C)
plot(ω, abs.(freqresp(G_yd, ω)[:]), xscale=:log, yscale=:log,
    title="Transfer function from load disturbances to output.")
