function simulate(sys::DelayLtiSystem, tspan; u=[], x0=Float64[], alg=MethodOfSteps(Tsit5()), kwargs...)
    ControlSystems.PartionedStateSpace(sys.P, 1, 1)

    nu = ControlSystems.ninputs(sys)
    nx = ControlSystems.nstates(sys)

    u = (u == []) ? t -> fill(0.0, nu) : u
    x0 = (x0 == []) ? fill(0.0, nx) : x0

    # FIXME: Only works for one single time delay
    # FIXME: Should error for non-zero D22 terms
    dde = function (dx, x, h, p, t)
        d_delayed = P.C2*h(p,t-sys.Tau[1]) + P.D21*u(t-sys.Tau[1])# + P.D22*(t-Tau[1]
        dx .= P.A*x  + P.B2*d_delayed + P.B1*u(t)
    end

    h_initial = (p, t) -> zeros(N)


    saved_values_y = SavedValues(Float64, Vector{Float64})
    cb = SavingCallback((x,t,integrator)->(P.C1*x + P.D11*u(t)), saved_values_y)

    prob = DDEProblem(dde, x0, h_initial, tspan, constant_lags=sys.Tau)

    sol = solve(prob, alg; saveat=0:0.02:tspan[2], callback=cb, kwargs...)

    t, x = sol.t, sol.u


    t, x
end

sys = feedback(1.0, ss(-1, 1, 1, 0) * delay(1.0))

@time t, x = simulate(sys, (0.0,4.0); u=t->[t>0 ? 1.0 : 0.0], saveat=0:0.02:5)

x1 = [x[k][1] for k=1:length(x)]
del = [zeros(50); x1[1:end-50]] # FIXME: Not the real delayed signal...

# The following does not work in general... just this specific problem instance
plot(t, ones(size(t)) .- x1)
