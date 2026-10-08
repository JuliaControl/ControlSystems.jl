using Plots

tf2 = tf([1/5,2],[1,1,1])
rts, Z, K = rlocus(tf2)
rts2, Z2, _ = rlocus(tf2; K=K)
@test size(rts) == size(rts2)

P, Q = numpoly(tf2)[], denpoly(tf2)[]
for (k, rs, rs2) = zip(eachrow(K), eachrow(rts), eachrow(rts2))
    # test that entries are solutions
    for r in rs
        @test isapprox((k[1]*P+Q)(r), 0, atol=1e-10)
    end
    @test isapprox(rs, rs2)
end
rlocusplot(tf2)
plot(rlocus(tf2))

# https://github.com/JuliaControl/ControlSystems.jl/issues/740
Nroots = -1.0 .+  [-1.732050807568877im, 1.732050807568877im]
Droots =  [0.0, -4.0, -6.0, -0.7 - 0.7141428428542851im, -0.7 + 0.7141428428542851im]
G = zpk(Nroots, Droots, 1.0)
rlocusplot(G, 200)

# The number of steps must not grow linearly with the distance the poles travel
sys = ss([-0.5 4000; -66.66666666666667 -875], [0.0; 6000.0;;], [9.549296585513721 0.0], 0)
r = rlocus(sys)
@test length(r.K) < 2000
P, Q = numpoly(tf(sys))[], denpoly(tf(sys))[]
for (k, rs) = zip(r.K, eachrow(r.roots))
    for p in rs
        @test abs((k*P+Q)(p)) <= 1e-8*abs(k*P(p))+1e-8*abs(Q(p))
    end
end

# The step-size control is independent of the time unit of the system
sys_ms = ss(sys.A/1000, sys.B/1000, sys.C, 0) # Time unit ms instead of s
r_ms = rlocus(sys_ms)
@test length(r_ms.K) == length(r.K)
@test r_ms.roots ≈ r.roots/1000

# Improper transfer function, the closed-loop system has more poles than the open-loop system
G = tf([1, 2, 3], [1, 1])
r = rlocus(G, 100)
@test size(r.roots, 2) == 2
@test 0 < r.K[1] < r.K[end] == 100
P, Q = numpoly(G)[], denpoly(G)[]
for (k, rs) = zip(r.K, eachrow(r.roots))
    for p in rs
        @test isapprox((k*P+Q)(p), 0, atol=1e-8*(abs(k*P(p))+abs(Q(p))))
    end
end
@test minimum(abs.(r.roots[1, :] .- poles(G)')) <= 1.01e-2 # One pole starts close to the open-loop pole
@test maximum(abs, r.roots[1, :]) > 10*maximum(abs, [poles(G); tzeros(G)]) # The other pole starts at a large magnitude
@test sort(r.roots[end, :], by=imag) ≈ sort(tzeros(G), by=imag) rtol=0.1 # All poles approach the zeros
rts, _ = ControlSystemsBase.getpoles(G, [0.0, 1.0, 2.0])
@test size(rts) == (3, 2)
plot(r)
