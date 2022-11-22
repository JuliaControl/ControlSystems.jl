using Plots

tf2 = tf([1/5,2],[1,1,1])
rts, Z, K = rlocus(tf2)
rts2, Z2, _ = rlocus(tf2; K=K)
P, Q = numpoly(tf2)[], denpoly(tf2)[]

@test size(rts) == size(rts2)
for (k, rs, rs2) = zip(eachrow(K), eachrow(rts), eachrow(rts2))
    # test that entries are solutions
    for r in rs
        @test isapprox((k[1]*P+Q)(r), 0, atol=1e-12)
    end
    @test isapprox(rs, rs2)
end

# https://github.com/JuliaControl/ControlSystems.jl/issues/740
Nroots = -1.0 .+  [-1.732050807568877im, 1.732050807568877im]
Droots =  [0.0, -4.0, -6.0, -0.7 - 0.7141428428542851im, -0.7 + 0.7141428428542851im]
G = zpk(Nroots, Droots, 1.0)
rlocusplot(G, 200)
