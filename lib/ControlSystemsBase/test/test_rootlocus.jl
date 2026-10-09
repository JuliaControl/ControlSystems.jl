using Plots

tf2 = tf([1/5,2],[1,1,1])
rts, Z, K = rlocus(tf2)
rts2, Z2, _ = rlocus(tf2; K=K)
@test size(rts) == size(rts2)

P, Q = numpoly(tf2)[], denpoly(tf2)[]
for (k, rs, rs2) = zip(eachrow(K), eachrow(rts), eachrow(rts2))
    # test that entries are solutions
    for r in rs
        @test isapprox((k[1]*P+Q)(r), 0, atol=1e-10*(1 + k[1]))
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

# Default gain: the poles that tend to infinity have a magnitude of at least ten times the largest magnitude among the open-loop poles and zeros, the remaining poles are close to the zeros
using ControlSystemsBase: default_rlocus_gain, rlocus_limits
for G in (tf([1/5,2],[1,1,1]), tf(1, [1, 3, 2, 0]), zpk([-1000], [-1, -10, -10], 1), tf([-1, 1], [1, 1, 1]), tf(sys))
    local r = rlocus(G)
    @test r.K[end] == default_rlocus_gain(G)
    ω0 = maximum(abs, [poles(G); tzeros(G)])
    z = tzeros(G)
    far = abs.(r.roots[end, :]) .> 10ω0
    @test count(far) == length(poles(G)) - length(z)
    @test all(minimum(abs.(p .- z)) <= 0.0101*max(abs(z[argmin(abs.(p .- z))]), 0.01ω0) for p in r.roots[end, .!far])
end
@test default_rlocus_gain(10tf(sys)) ≈ default_rlocus_gain(tf(sys))/10
@test default_rlocus_gain(sys_ms) ≈ default_rlocus_gain(sys)
@test_throws ErrorException rlocus(ssrand(2, 2, 3))

# Default plot limits contain the open-loop poles and the origin, and are not determined by the poles at the final gain
r = rlocus(sys)
r_ms = rlocus(sys_ms)
xl, yl = rlocus_limits(r)
ω0 = maximum(abs, poles(sys))
@test all(xl[1] < real(p) < xl[2] && yl[1] < imag(p) < yl[2] for p in [poles(sys); 0])
@test xl[2] - xl[1] < 3ω0 && yl[2] - yl[1] < 6ω0
xl_ms, yl_ms = rlocus_limits(r_ms)
@test all(xl_ms .≈ xl ./ 1000) && all(yl_ms .≈ yl ./ 1000)

# Discrete-time system, the limits contain the unit circle
Gd = zpk([-0.5], [1, 0.6], 0.1, 0.1)
rd = rlocus(Gd)
xl, yl = rlocus_limits(rd)
@test xl[1] < -1 && xl[2] > 1 && yl[1] < -1 && yl[2] > 1
plot(rd)
rlocusplot(Gd)
