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
