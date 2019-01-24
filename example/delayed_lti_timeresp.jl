using Plots
plotly()

sys = feedback(1.0, ss(-1.0, 2, 1, 0) * (delay(2.0) + delay(3.0) + delay(2.5)))

t = 0:0.02:8
@time t, x, y = lsim(sys, t; u=t->[t>=0 ? 1.0 : 0.0])

plot(t, y')
gui()
