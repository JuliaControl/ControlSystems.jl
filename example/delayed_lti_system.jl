using Plots

# Frequency domain analysis of a simple system with time delays

s = tf("s")

P = delay(0.2) * ss(1/(s+1))

K = 3; Ti = 0.3;
C = ss(K*(1 + 1/s))

ω = exp10.(LinRange(-2,2,500))

nyquistplot(P*C, ω, xlims=(-2,1), ylims=(-2,2), gaincircles=false)

G_yd = feedback(P, C)
bodeplot(G_yd, ω, plotphase=false, title="Transfer function from load disturbances to output.")
