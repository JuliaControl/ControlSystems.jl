@testset "test_synthesis" begin
P = tf(1.,[1.,1])
C = tf([1.,1],[1.,0])
L = P*C
Lsys = ss(L)


B = [1]
A = [1,1]
R = [1,1]
S = [1]
T = [1]
@testset "minreal + feedback" begin
@test isapprox(minreal(feedback(P,C),1e-5), tf([1,0],[1,2,1]), rtol = 1e-5)
@test isapprox(numpoly(minreal(feedback(L),1e-5))[1].coeffs, numpoly(tf(1,[1,1]))[1].coeffs)# This test is ugly, but numerical stability is poor for minreal
@test feedback2dof(B,A,R,S,T) == tf(B.*T, conv(A,R) + [0;0;conv(B,S)])
@test feedback2dof(P,R,S,T) == tf(B.*T, conv(A,R) + [0;0;conv(B,S)])
@test isapprox(poles(minreal(tf(feedback(Lsys)),1e-5)) , poles(minreal(feedback(L),1e-5)), atol=1e-5)

Pint = tf(1,[1,1])
Cint = tf([1,1],[1,0])
Lint = P*C

@test isapprox(minreal(feedback(Pint,Cint),1e-5), tf([1,0],[1,2,1]), rtol = 1e-5) # TODO consider keeping minreal of Int system Int
@test isapprox(numpoly(minreal(feedback(Lint),1e-5))[1].coeffs, numpoly(tf(1,[1,1]))[1].coeffs)# This test is ugly, but numerical stability is poor for minreal
@test isapprox(poles(minreal(tf(feedback(Lsys)),1e-5)) , poles(minreal(feedback(L),1e-5)), atol=1e-5)

@test feedback(ss(1),ss(1)) == ss(0.5)
@test_throws ErrorException feedback(ss([1 0; 0 1], ones(2,2), ones(1,2),0))

# Test Feedback Issue: 163
g = tf([1],[1,1])
gfb = feedback(g)
gfb2 = tf(feedback(ss(g)))
@test hinfnorm(gfb - gfb2)[1] <= 1e-14

# Test more feedback
s = tf("s")
ftf = 1.0*(2s+3)/((5s+7)*(11s+13))
fzpk = zpk(ftf)

ffb = feedback(fzpk)            # Zpk feedback
ffb2 = zpk(feedback(ss(fzpk)))  # ss feedback + ss conversion
ffb3 = feedback(ftf)            # tf feedback
ffb4 = feedback(ss(ftf))        # ss feedback
ffb5 = minreal(fzpk/(1+fzpk))   # Zpk feedback manual

z1,p1,k1 = zpkdata(ffb)
z2,p2,k2 = zpkdata(ffb2)
z3,p3,k3 = zpkdata(ffb3)
z4,p4,k4 = zpkdata(ffb4)
z5,p5,k5 = zpkdata(ffb5)
@test sort(real.(z1[1])) ≈ sort(real.(z2[1])) ≈ sort(real.(z3[1])) ≈ sort(real.(z4[1])) ≈ sort(real.(z5[1]))
@test sort(real.(p1[1])) ≈ sort(real.(p2[1])) ≈ sort(real.(p3[1])) ≈ sort(real.(p4[1])) ≈ sort(real.(p5[1]))
@test k1 ≈ k2 ≈ k3 ≈ k4 ≈ k5
end



@testset "place" begin
    sys = ss(-4, 2, 3, 0)
    A, B, C, _ = ssdata(sys)

    @test place(A, B, [-10]) == [3][:,:]
    @test place(A, B, [-10], :c) == [3][:,:]
    @test place(A, C, [-10], :o) == [2][:,:]

    A = [0 1; 0 0]
    B = [0; 1]
    C = [1 0]
    sys = ss(A, B, C, 0)

    @test place(A, B, [-1.0, -1]) ≈ [1 2]
    @test place(sys, [-1.0, -1]) ≈ [1 2]
    @test place(A, B, [-1.0, -1], :c) ≈ [1 2]
    @test place(sys, [-1.0, -1], :c) ≈ [1 2]
    @test place(A, C, [-2.0, -2], :o) ≈ [4; 4]
    @test place(sys, [-2.0, -2], :o) ≈ [4; 4]

    @test place(A, B, [-2 + im, -2 - im]) ≈ [5 4]
    @test place(A, C, [-4 + 2im, -4 - 2im], :o) ≈ [8; 20]

    A = ones(3,3) - diagm([3, 4, 5])
    B = [1; 0; 2]
    C = [1 1 0]
    @test place(A, B, [-2 + 2im, -2 - 2im, -4]) ≈ [-2.6 5.2 0.8]
    @test place(A, C, [-2 + 3im, -2 - 3im, -4], :o) ≈ [11; -12; 1]
end


@testset "acker" begin
Random.seed!(0)
A = randn(3,3)
B = randn(3,1)
p = [3.0,2,1]
K = ControlSystemsBase.acker(A,B,p)
@test ControlSystemsBase.eigvalsnosort(A-B*K) ≈ p

p = [-1+im, -1-im, -1]
K = ControlSystemsBase.acker(A,B,p)
@test ControlSystemsBase.eigvalsnosort(A-B*K) ≈ p
end


@testset "MIMO place" begin

    function allin(a, b; tol = 1e-2)
        all(minimum(abs, a .- b', dims=2)[:] .< tol)
    end

    for i = 1:10
        # B smaller than A
        sys = ssrand(2,2,3)
        # @show cond(gram(sys, :c))
        (; A, B) = sys
        p = [-1.0, -2, -3]
        L = place(A, B, p, verbose=false)
        @test eltype(L) <: Real
        @test allin(eigvals(A - B*L), p)

        # p = [-3.0, -1-im, -1+im]
        # L = place(A, B, p)
        # @test eltype(L) <: Real

        # cond(gram(sys, :c))
        # @test allin(eigvals(A - B*L), p)


        # B same size as A
        sys = ssrand(2,3,3)
        (; A, B) = sys
        p = [-1.0, -2, -3]
        L = place(A, B, p)
        @test eltype(L) <: Real
        @test allin(eigvals(A - B*L), p)

        p = [-3.0, -1-im, -1+im]
        L = place(A, B, p)
        @test eltype(L) <: Real
        @test allin(eigvals(A - B*L), p)

        # deadbeat
        A = [0 1; 0 0]
        B = I(2)
        sys = ss(A, B, I, 0)
        sysd = c2d(sys, 0.1)

        p = [0,0]
        L = place(sysd, p)
        @test eltype(L) <: Real
        @test allin(eigvals(sysd.A - sysd.B*L), p)

        # B of size 1
        sys = ssrand(1,1,3)
        (; A, B) = sys
        p = [-1.0, -2, -3]
        L = ControlSystemsBase.place_knvd(A, B, p)
        @test eltype(L) <: Real
        @test allin(eigvals(A - B*L), p)

    end

    A = [
        -0.1094 0.0628 0 0 0
        1.306 -2.132 0.9807 0 0
        0 1.595 -3.149 1.547 0
        0 0.0355 2.632 -4.257 1.855
        0 0.00227 0 0.1636 -0.1625
    ]
    B = [
        0 0.0638 0.0838 0.1004 0.0063
        0 0 -0.1396 -0.206 -0.0128
    ]'
    # p = [-0.07732, -0.01423, -0.8953, -2.841, -5.982]
    p = [-0.2, -0.5, -1, -1+im, -1-im]
    L = place(A, B, p; verbose=true) # with verbose, should prind cond(X) ≈ 39.4 which it does 
    @test allin(eigvals(A - B*L), p; tol=0.025) # Tolerance from paper
    # norm(L)

    ## Rank deficient multi-input B
    P = let
        tempA = [0.0 1.0; -4.0 -1.2]
        tempB = [0.0 0.0; 4.0 4.0]
        tempC = [1.0 0.0]
        tempD = [0.0 0.0]
        ss(tempA, tempB, tempC, tempD)
    end
    F = place(P, [-2, -2], verbose=true)
    @test allin(eigvals(P.A - P.B*F), [-2, -2])


    P = let
        tempA = [0.0 1.0; -4.0 -1.2]
        tempB = randn(2, 1) * randn(1, 10) # Rank deficient
        tempC = [1.0 0.0]
        tempD = zeros(1, 10)
        ss(tempA, tempB, tempC, tempD)
    end
    F = place(P, [-2, -2], verbose=true)
    @test allin(eigvals(P.A - P.B*F), [-2, -2])

end

@testset "LQR" begin
    Ts = 0.1
    A = [1 Ts; 0 1]
    B = [0;1] # Note B is vector, B'B is scalar, but compatible with I
    C = [1 0]
    Q = I
    R = I
    L = lqr(Discrete, A,B,Q,R)
    @test L ≈ [0.5890881713787511 0.7118839434795103]
    L, S, p = lqr(Discrete, A,B,Q,R, extra=Val(true))
    @test p ≈ eigvals(A-B*L)
    sys = ss(A,B,C,0,Ts)
    L = lqr(sys, Q, R)
    @test L ≈ [0.5890881713787511 0.7118839434795103]

    L = lqr(sys, Q, R)
    @test L ≈ [0.5890881713787511 0.7118839434795103]

    B = reshape(B,2,1)  # Note B is matrix, B'B is compatible with I
    L = lqr(Discrete, A,B,Q,R)
    @test L ≈ [0.5890881713787511 0.7118839434795103]

    Q = eye_(2)
    R = eye_(1)
    L = lqr(Discrete, A,B,Q,R)
    @test L ≈ [0.5890881713787511 0.7118839434795103]

    B = [0;1]   # Note B is vector, B'B is scalar 
    Q = eye_(2)
    R = eye_(1)
    L ≈ lqr(Discrete, A,B,Q,R)
    #L ≈ [0.5890881713787511 0.7118839434795103]
end

@testset "issue #304" begin
    @test feedback(ss(0.5, 1, 1, 0, 1), 1) == ss(-0.5, 1, 1, 0, 1)
    @test feedback(1, ss(0.5, 1, 1, 0, 1)) == ss(-0.5, 1, -1, 1, 1)
    @test [ss(1,1,1,0,1) 1] == ss(1, [1 0], 1, [0 1], 1)
end

end

@testset "alpha_beta / alpha_beta_gamma" begin
    Ts = 0.1
    kalata(a) = 2 * (2 - a) - 4 * sqrt(1 - a)

    @testset "shape and defaults" begin
        s2 = alpha_beta(0.5, Ts)
        s3 = alpha_beta_gamma(0.5, Ts)
        @test size(s2) == (2, 1)
        @test size(s3) == (3, 1)
        @test s2.Ts == Ts && s3.Ts == Ts
        @test iszero(s2.D) && iszero(s3.D)          # the state is the estimate
        @test s2.C == I && s3.C == I

        # The documented default gains.
        @test kalata(0.5) ≈ 0.1715728752538097
        @test kalata(0.5)^2 / (2 * 0.5) ≈ 0.029437251522859908
        # ... and they are what the systems are actually built with: B == K.
        @test s2.B ≈ [0.5; kalata(0.5) / Ts;;]
        @test s3.B ≈ [0.5; kalata(0.5) / Ts; kalata(0.5)^2 / (2 * 0.5) / Ts^2;;]
    end

    @testset "reproduces the documented recurrence" begin
        a, b, g = 0.4, 0.2, 0.05
        t = 0:Ts:2
        us = collect(float.(t))                      # a unit-slope ramp

        x = 0.0; v = 0.0; X = Float64[]; V = Float64[]
        for u in us
            pred = x + Ts * v
            r = u - pred
            x = pred + a * r
            v = v + (b / Ts) * r
            push!(X, x); push!(V, v)
        end
        y = lsim(alpha_beta(a, Ts; beta = b), reshape(us, 1, :), t).y
        # The estimate that absorbed us[k] appears at output index k+1.
        @test y[1, 2:end] ≈ X[1:end-1]
        @test y[2, 2:end] ≈ V[1:end-1]

        x = 0.0; v = 0.0; ac = 0.0; X = Float64[]; V = Float64[]; A = Float64[]
        for u in us
            pred = x + Ts * v + Ts^2 / 2 * ac
            predv = v + Ts * ac
            r = u - pred
            x = pred + a * r
            v = predv + (b / Ts) * r
            ac = ac + (g / Ts^2) * r
            push!(X, x); push!(V, v); push!(A, ac)
        end
        y = lsim(alpha_beta_gamma(a, Ts; beta = b, gamma = g), reshape(us, 1, :), t).y
        @test y[1, 2:end] ≈ X[1:end-1]
        @test y[2, 2:end] ≈ V[1:end-1]
        @test y[3, 2:end] ≈ A[1:end-1]
    end

    @testset "tracks a ramp" begin
        t = 0:Ts:5
        res = lsim(alpha_beta(0.5, Ts), (x, t) -> [t], t)
        @test res.y[2, end] ≈ 1 atol = 1e-6          # rate converges to the slope
        res = lsim(alpha_beta_gamma(0.5, Ts), (x, t) -> [t], t)
        @test res.y[2, end] ≈ 1 atol = 1e-3          # third order settles slower on a ramp
        @test res.y[3, end] ≈ 0 atol = 1e-3          # and the acceleration to zero
    end

    @testset "equals observer_filter of the integrator chain" begin
        for (n, f) in ((2, alpha_beta), (3, alpha_beta_gamma))
            A = n == 2 ? [1 Ts; 0 1.0] : [1 Ts Ts^2/2; 0 1 Ts; 0 0 1.0]
            B = n == 2 ? [Ts^2/2; Ts] : [Ts^3/6; Ts^2/2; Ts]
            C = n == 2 ? [1.0 0] : [1.0 0 0]
            sys = ss(A, B, C, 0, Ts)
            a, b, g = 0.4, 0.2, 0.05
            K = n == 2 ? [a; b/Ts;;] : [a; b/Ts; g/Ts^2;;]
            filt = n == 2 ? f(a, Ts; beta = b) : f(a, Ts; beta = b, gamma = g)
            ref = observer_filter(sys, K; output_state = true)
            @test filt.A ≈ ref.A
            @test filt.B ≈ ref.B[:, 2:2]             # observer_filter also takes u; this has only y
        end
    end

    @testset "critical damping" begin
        # The recipes in the extended help place every error pole on the real axis at s. The
        # eigenvalue is defective, so it is only conditioned to about sqrt(eps).
        for s in (0.95, 0.9, 0.8, 0.5, 0.2)
            p = eigvals(alpha_beta(1 - s^2, Ts; beta = (1 - s)^2).A)
            @test all(z -> isapprox(z, s; atol = 1e-6), p)

            p = eigvals(alpha_beta_gamma(1 - s^3, Ts;
                                         beta = 1.5 * (1 - s)^2 * (1 + s), gamma = (1 - s)^3).A)
            @test all(z -> isapprox(z, s; atol = 1e-4), p)
        end

        # Independent of the sample rate: these gains are dimensionless in this parameterization.
        for Ts2 in (1.0, 0.01, 1e-4)
            p = eigvals(alpha_beta(1 - 0.9^2, Ts2; beta = (1 - 0.9)^2).A)
            @test all(z -> isapprox(z, 0.9; atol = 1e-6), p)
        end

        # The defaults are a different tuning and always leave a complex pair, which is why the
        # recipes have to set every gain rather than just `alpha`.
        for a in (0.2, 0.5, 0.9)
            @test any(z -> abs(imag(z)) > 1e-6, eigvals(alpha_beta(a, Ts).A))
            @test any(z -> abs(imag(z)) > 1e-6, eigvals(alpha_beta_gamma(a, Ts).A))
        end
        # The numbers the docstring quotes for s = 0.9.
        @test all((1 - 0.9^3, 1.5 * 0.1^2 * 1.9, 0.1^3) .≈ (0.271, 0.0285, 0.001))
        @test kalata(0.271) ≈ 0.0427 atol = 1e-4
        @test kalata(0.271)^2 / (2 * 0.271) ≈ 0.0034 atol = 1e-4
    end

    @testset "argument checking" begin
        @test_throws ArgumentError alpha_beta(0.0, Ts)
        @test_throws ArgumentError alpha_beta(1.0, Ts)
        @test_throws ArgumentError alpha_beta(0.5, 0.0)
        @test_throws ArgumentError alpha_beta_gamma(-0.1, Ts)
        @test_throws ArgumentError alpha_beta_gamma(0.5, -1.0)
    end
end
