@testset "Transfer-function minreal tolerances" begin
    # Inspect the stored roots: poles(::TransferFunction) can itself remove
    # common factors through MatrixPencils and would hide a missed cancellation.
    function roots_of(G)
        z, p, _ = zpkdata(G)
        return only(z), only(p)
    end

    @testset "Slow dynamics (#1081)" begin
        G = tf([1, 1.5e-8], [1, 1e-8])
        R = minreal(G)
        @test dcgain(R)[] ≈ 1.5
        @test length(last(roots_of(R))) == 1
        @test freqresp(R, [0, 1e-8, 1]) ≈ freqresp(G, [0, 1e-8, 1])
        for a in (1e-12, 1, 1e12)
            H = minreal(time_scale(G, a))
            @test length(last(roots_of(H))) == 1
            @test dcgain(H)[] ≈ 1.5
            @test freqresp(H, a.*[0, 1e-8, 1]) ≈ freqresp(G, [0, 1e-8, 1])
        end
        @test dcgain(zpk([], [-1e-8], 1e-8) + 1)[] ≈ 2

        # An unrelated fast pole must not enlarge the cancellation threshold.
        G = zpk([-1.5e-8], [-1e-8, -1e12], 1e12)
        @test roots_of(minreal(G)) == roots_of(G)
        @test dcgain(minreal(G))[] ≈ 1.5

        # Explicit tolerances retain their absolute meaning, including zero.
        @test isempty(last(roots_of(minreal(zpk([-1.5e-8], [-1e-8], 1), 1e-8))))
        G = zpk([-1.0], [-1.0], 1)
        @test roots_of(minreal(G, 0)) == roots_of(G)
    end

    @testset "Root types and scales" begin
        for T in (Float32, Float64, BigFloat), scale in T.((1e-12, 1, 1e12))
            δ = eps(T)*4
            G = zpk(T[-scale*(1+δ)], T[-scale], one(T))
            @test isempty(last(roots_of(@inferred minreal(G))))
            @test !isempty(last(roots_of(minreal(G, zero(T)))))
            G = zpk(T[-1.5scale], T[-scale], one(T))
            @test roots_of(minreal(G)) == roots_of(G)
            G = tf(T[1, 1.5scale], T[1, scale])
            @test dcgain(minreal(G))[] ≈ T(1.5)
            @test length(last(roots_of(minreal(G)))) == 1

            # An exact common root, including repeated roots and the origin.
            G = zpk(T[0, -scale, -scale], T[0, -scale, -scale, -2scale], one(T))
            @test roots_of(minreal(G)) == (T[], T[-2scale])
            G = zpk(T[0], T[-scale], one(T))
            @test roots_of(minreal(G)) == roots_of(G)
            G = zpk(T[scale], T[-scale], one(T))
            @test roots_of(minreal(G)) == roots_of(G)
        end

        # Distances must not lose the type's precision to Float64 defaults.
        setprecision(256) do
            G = zpk(BigFloat[-1- big"1e-30"], BigFloat[-1], big"1")
            @test roots_of(minreal(G)) == roots_of(G)
        end
        for T in (Float32, Float64)
            scale = nextfloat(zero(T))
            G = zpk(T[-2scale], T[-scale], one(T))
            @test roots_of(minreal(G)) == roots_of(G)
            # abs(complex root) overflows although its components are finite.
            p = complex(floatmax(T), floatmax(T))
            G = zpk([p/2], [p], complex(one(T)))
            @test roots_of(minreal(G)) == roots_of(G)
            G = zpk([p*(1-4eps(T))], [p], complex(one(T)))
            @test isempty(last(roots_of(minreal(G))))
        end
    end

    @testset "Conjugates, ordering, and time evolution" begin
        for T in (Float32, Float64), scale in T.((1e-12, 1, 1e12))
            c = complex(-scale, scale)
            δ = eps(T)*4
            G = zpk([c*(1+δ), conj(c)*(1+δ), -2scale],
                    [c, conj(c), -3scale], one(T))
            R = minreal(G)
            @test roots_of(R) == (Complex{T}[-2scale], Complex{T}[-3scale])
            # Evaluate in Float64 to avoid overflow in Float32 cubic products.
            w = Float64.(scale.*T[0, 0.5, 1, 2, 10])
            @test freqresp(R, w) ≈ freqresp(G, w) rtol=sqrt(eps(T))
            # Retain distinct conjugate pairs in an actual real polynomial TF.
            G = tf(T[1, 2scale, 2scale^2], T[1, 2scale, 3scale^2])
            @test length(last(roots_of(minreal(G)))) == 2

            # Multiple nearby roots must consume a matching zero only once.
            G = zpk(T[-scale, -scale*(1+8sqrt(eps(T)))], T[-scale, -scale], one(T))
            R = minreal(G)
            @test length(first(roots_of(R))) == length(last(roots_of(R))) == 1
        end
        for Ts in (0.01, 1, 100)
            G = zpk([0.25, 0.5], [0.25, 0.75], 2.0, Ts)
            R = minreal(G)
            @test R.Ts == Ts
            @test roots_of(R) == ([0.5], [0.75])
            @test freqresp(R, [0, 0.1, 1]./Ts) ≈ freqresp(G, [0, 0.1, 1]./Ts)
        end
        f = zpk([-2e6], [-1e6, -3e6, -4e6], 1e12)
        @test length(last(roots_of(f+f))) == 3
        @test evalfr(f+f, 1e6im) ≈ 2evalfr(f, 1e6im)

        # Apply the policy independently to every matrix entry.
        G = [tf([1, 1.5e-8], [1, 1e-8]) tf([1, 2, 1], [1, 3, 2])]
        _, ps, _ = zpkdata(minreal(G))
        @test length.(ps) == [1 1]
        @test dcgain(minreal(G)) ≈ dcgain(G)
    end

    @testset "Independent polynomial common factors" begin
        rng = MersenneTwister(1081)
        for T in (Float32, Float64), scale in T.((1e-12, 1, 1e12)), _ in 1:20
            # Separated simple roots keep conditioning out of this regression.
            c, a, b = T(0.5)+rand(rng,T), T(3)+rand(rng,T), T(6)+rand(rng,T)
            G = tf(T[1, (c+a)*scale, c*a*scale^2],
                   T[1, (c+b)*scale, c*b*scale^2])
            R = minreal(G)
            @test length(first(roots_of(R))) == length(last(roots_of(R))) == 1
            @test dcgain(R)[] ≈ a/b rtol=sqrt(eps(T))
            @test freqresp(R, scale.*T[0, 1, 10]) ≈
                  freqresp(G, scale.*T[0, 1, 10]) rtol=sqrt(eps(T))
            # Independently evaluate the original rounded coefficients, without
            # computing their roots or using the same frequency-response path.
            s = Complex{BigFloat}(0, BigFloat(scale))
            reference = evalpoly(s, reverse(BigFloat.(only(numvec(G))))) /
                        evalpoly(s, reverse(BigFloat.(only(denvec(G)))))
            @test evalfr(R, s)[] ≈ reference rtol=sqrt(eps(T))
        end
    end
end
