using Plots
plotly()

sys = feedback(1.0, ss(-1.0, 2, 1, 0) * (delay(2.0) + delay(3.0) + delay(2.5)))

t = 0:0.02:8
@time lsim(sys, t-> [t>=0 ? 1.0 : 0.0], t)
@time t, x, y = lsim(sys, [1.0], t)
@time t, x, y = lsim(sys, (out, t) -> (out .= (t>=0 ? 1.0 : 0.0)), t)
@time t, x, y = lsim(sys, (out, t) -> (out[1] = (t>=0 ? 1.0 : 0.0)), t)

function u0(out,t)
    if t > 0
        out[1] = 1
    else
        out[1] = 0
    end
    return
end

@time t, x, y = lsim(sys, t; u=u0)

plot(t, y')
gui()

@time ControlSystems._lsim(sys, t, (out,t)-> out[1] = (t>=0 ? 1.0 : 0.0))
