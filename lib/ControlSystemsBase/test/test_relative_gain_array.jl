@testset "Relative gain array numerics" begin
    rng = MersenneTwister(1081)

    # Independent UC reference: double-center log magnitudes in one step
    # (full support), then use a known rank-two factorization in 256 bits.
    A0 = [1 8 17; -4 6 8; -8 1 -6] # column 3 = column 1 + 2column 2
    Rref = setprecision(256) do
        A = BigFloat.(A0)
        L = log.(abs.(A))
        u = -sum(L; dims=2) / 3
        v = -sum(L; dims=1) / 3 .+ sum(L) / 9
        S = exp.(u) .* A .* exp.(v)
        B = S[:, 1:2]
        C = BigFloat[1 0 exp(v[3]-v[1]); 0 1 2exp(v[3]-v[2])]
        SI = C' * inv(C*C') * inv(B'*B) * B'
        S .* transpose(SI)
    end
    for T in (Float32, Float64)
        τ = 100eps(T)
        A = T.(A0)
        @test eltype(@inferred relative_gain_array(A)) == T
        for a in T.((1, 1e-3, 1e-12, 1e12))
            @test relative_gain_array(a*A) ≈ Rref rtol=τ
        end
        # A rank-one matrix stays rank one, including with zero rows/columns.
        for (m,n) in ((2,3), (3,2), (3,3))
            A = T.(collect(1:m) * collect(1:n)')
            @test relative_gain_array(A) ≈ fill(T(1)/(m*n),m,n) rtol=τ
            Z = [zeros(T,1,n+1); zeros(T,m) A]
            @test relative_gain_array(Z) ≈ [zeros(T,1,n+1); zeros(T,m) fill(T(1)/(m*n),m,n)] rtol=τ
        end
        @test relative_gain_array(zeros(T,2,3)) == zeros(T,2,3)
        @test relative_gain_array(zeros(T,0,3)) == zeros(T,0,3)
        @test relative_gain_array(zeros(T,3,0)) == zeros(T,3,0)
        @test relative_gain_array(T[1 2 0 0; 2 4 0 0; 0 0 1 3]) ≈
              T[0.25 0.25 0 0; 0.25 0.25 0 0; 0 0 0.5 0.5] rtol=τ

        # Full-rank square, tall, wide and rank-deficient matrices, with
        # independent signed input/output scalings and permutations.
        for A in (T[4 1 2; 1 5 1; 2 1 6], T[1 2 0; 0 1 1],
                  T[1 0; 2 1; 0 1], T.(A0))
            m,n = size(A)
            R = relative_gain_array(A)
            d = T.((-1).^(1:m) .* exp10.(range(-8,8;length=m)))
            e = T.((-1).^(1:n) .* exp10.(range(8,-8;length=n)))
            for B in (T(1e-12)*A, d .* A, A .* e', d .* A .* e')
                @test relative_gain_array(B) ≈ R rtol=τ
            end
            p = randperm(rng,m); q = randperm(rng,n)
            @test relative_gain_array(A[p,q]) ≈ R[p,q] rtol=τ
            @test relative_gain_array(transpose(A)) ≈ transpose(R) rtol=τ
            @test sum(R) ≈ rank(A) rtol=τ
        end

        # Random matrices checked against the inverse of the supplied input
        # values, not against their pre-rounding generating matrix.
        for _ in 1:10
            A = randn(rng,T,4,4) + 4I
            R = relative_gain_array(A)
            B = BigFloat.(A)
            @test R ≈ B .* transpose(inv(B)) rtol=τ
            d = T.(exp10.(range(-6,6;length=4)))
            e = reverse(d)
            @test relative_gain_array(d .* A .* e') ≈ R rtol=τ
            @test vec(sum(R;dims=1)) ≈ ones(T,4) rtol=τ
            @test vec(sum(R;dims=2)) ≈ ones(T,4) rtol=τ
        end

        # A resolvable independent direction must not be discarded merely
        # because it is below sqrt(eps). The exact 2x2 RGA is available.
        δ = T === Float32 ? T(2)^(-12) : T(2)^(-40)
        A = T[1 1; 1 1+δ]
        B = BigFloat.(A)
        @test relative_gain_array(A) ≈ B .* transpose(inv(B)) rtol=τ
        @test maximum(abs,relative_gain_array(A)) > inv(2δ)

        CT = Complex{T}
        for A in (CT[1+im 2-im; 3im 4], CT[1+im 2 3im; 2+2im 4 6im])
            R = relative_gain_array(A)
            @test eltype(@inferred relative_gain_array(A)) == CT
            if size(A,1) == size(A,2)
                B = Complex{BigFloat}.(A)
                @test R ≈ B .* transpose(inv(B)) rtol=τ
            else
                @test R ≈ fill(T(1)/6,size(A)) rtol=τ
            end
            d = CT[im, -2]
            e = CT.((1:size(A,2)) .* (1+im))
            @test relative_gain_array(d .* A .* transpose(e)) ≈ R rtol=τ
            @test relative_gain_array(conj.(A)) ≈ conj.(R) rtol=τ
        end
    end

    # Range-safe equilibration must preserve a correct RGA, not just
    # invariance of an incorrect result. The raw inverse loses an entry here.
    A = [1e300 1e-300; 1e300 2e-300]
    @test relative_gain_array(A) ≈ [2 -1; -1 2]
    @test relative_gain_array(Matrix(Diagonal([1e300,1e-300]))) ≈ Matrix(1.0I,2,2)
    @test relative_gain_array(fill(nextfloat(0.0),2,3)) ≈ fill(1/6,2,3)
    @test_throws ErrorException relative_gain_array([1e300 1e300 1e-300; 1e300 1e300 1e-300; 1e-300 1e-300 1e300])
    A = ComplexF64[floatmax()*(1+im) floatmax(); floatmax()*(1+im) floatmax()]
    @test relative_gain_array(A) ≈ fill(0.25+0im,2,2)
    @test relative_gain_array([1 2; 3 5]) ≈ [-5 6; 6 -5]
    @test relative_gain_array(sparse([1.0 2 0; 0 1 1])) ≈ relative_gain_array([1.0 2 0; 0 1 1])
    @test relative_gain_array(BigFloat[1 2; 3 5]) ≈ BigFloat[-5 6; 6 -5]
    @test_throws MethodError relative_gain_array(BigFloat[1 2 3; 2 4 6])

    # Sparse support gives slow alternating projections. A finite cap must
    # fail explicitly instead of returning an unconverged UC result.
    A = [2.0 3 0; 0 4 5]
    @test_throws ErrorException relative_gain_array(A; maxiter=1)
    @test all(isfinite,relative_gain_array(A))
    @test relative_gain_array(1e100A) ≈ relative_gain_array(A)
    A = zeros(32,33)
    for i in 1:32
        A[i,i] = 2
        A[i,i+1] = 3
    end
    @test_throws ErrorException relative_gain_array(A)
    S = Float64.(A .!= 0) # On this tree support the UC magnitudes are all one.
    @test relative_gain_array(A; maxiter=50000) ≈ S .* transpose(pinv(S))
    @test_throws ArgumentError relative_gain_array(A; tol=0)
    @test_throws ArgumentError relative_gain_array(A; tol=-1)
    @test_throws ArgumentError relative_gain_array(A; tol=NaN)
    @test_throws ArgumentError relative_gain_array(A; maxiter=0)
    @test_throws ArgumentError relative_gain_array([1.0 NaN])
    @test_throws ArgumentError relative_gain_array([1.0 Inf])
end
