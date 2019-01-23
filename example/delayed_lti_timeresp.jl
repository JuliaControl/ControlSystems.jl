using DifferentialEquations






#sys_d = c2d(sys.P)

#function sim_discr(sys::DelayLtiSystem)

sys = feedback(2.0, ss(-1.0, 2, 1, 0) * delay(3.0))

@time t, x, saved_y = ControlSystems.simulate2(sys, (0.0,8.0); u=t->[t>0 ? 1.0 : 0.0], saveat=0:0.02:5)


saved_y
x1 = [x[k][1] for k=1:length(x)]


del = [zeros(50); x1[1:end-50]] # FIXME: Not the real delayed signal...

# The following does not work in general... just this specific problem instance
plot(t, ones(size(t)) .- x1)
