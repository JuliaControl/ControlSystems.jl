using PrecompileTools


PrecompileTools.@setup_workload begin
    P = ssrand(1,1,2, proper=true)
    PrecompileTools.@compile_workload begin
        Pd = c2d(P, 0.1)
        for P in (P, Pd)
            bode(P)
            # nyquist(P)
            step(P)
            feedback(P)
            feedback(P, P)
            minreal(P)
            balance_statespace(P)
            P*P
            P+P
            2.0*P
            # hinfnorm(P)
            poles(P)
            tzeros(P)
        end
        G = tf(1.0, [1.0, 1])
        ss(G)
        ss(G, minimal=true)
        G = tf(1.0, [1.0, 1], 1)
        ss(G)
        ss(G, minimal=true)

        # Pdel = P*delay(1.0)
        # pade(Pdel, 2)
    end
end