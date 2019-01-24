plotly()

sys = feedback(1.0, ss(-1.0, 2, 1, 0) * (delay(2.0) + delay(3.0)))

@time t, x, saved_y = ControlSystems.simulate5(sys, (0.0,8.0); u=t->[t>0 ? 1.0 : 0.0], saveat=0:0.02:8)

y = hcat([saved_y.saveval[k][1] for k=1:length(t)]...)
d = hcat([saved_y.saveval[k][2] for k=1:length(t)]...)

for k=1:length(sys.Tau)
    N_del = Integer(50*sys.Tau[k])
    dk = [zeros(N_del); d[k, 1:end-N_del]]

    for j=1:length(t)
        y[:, j] .+= sys.P.D12[:, k] * dk[j]
    end
end

plot(t, y')
gui()
