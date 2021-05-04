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
@test isapprox(pole(minreal(tf(feedback(Lsys)),1e-5)) , pole(minreal(feedback(L),1e-5)), atol=1e-5)

Pint = tf(1,[1,1])
Cint = tf([1,1],[1,0])
Lint = P*C

@test isapprox(minreal(feedback(Pint,Cint),1e-5), tf([1,0],[1,2,1]), rtol = 1e-5) # TODO consider keeping minreal of Int system Int
@test isapprox(numpoly(minreal(feedback(Lint),1e-5))[1].coeffs, numpoly(tf(1,[1,1]))[1].coeffs)# This test is ugly, but numerical stability is poor for minreal
@test isapprox(pole(minreal(tf(feedback(Lsys)),1e-5)) , pole(minreal(feedback(L),1e-5)), atol=1e-5)


@test_throws ErrorException feedback(ss(1),ss(1))
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


@testset "acker" begin
Random.seed!(0)
A = randn(3,3)
B = randn(3,1)
p = [3.0,2,1]
K = ControlSystems.acker(A,B,p)
@test ControlSystems.eigvalsnosort(A-B*K) ≈ p

p = [-1+im, -1-im, -1]
K = ControlSystems.acker(A,B,p)
@test ControlSystems.eigvalsnosort(A-B*K) ≈ p
end

@testset "LQR" begin
    Ts = 0.1
    A = [1 Ts; 0 1]
    B = [0;1] # Note B is vector, B'B is scalar, but compatible with I
    C = [1 0]
    Q = I
    R = I
    L = dlqr(A,B,Q,R)
    @test L ≈ [0.5890881713787511 0.7118839434795103]
    sys = ss(A,B,C,0,Ts)
    L = lqr(sys, Q, R)
    @test L ≈ [0.5890881713787511 0.7118839434795103]

    L = dlqr(sys, Q, R)
    @test L ≈ [0.5890881713787511 0.7118839434795103]

    B = reshape(B,2,1)  # Note B is matrix, B'B is compatible with I
    L = dlqr(A,B,Q,R)
    @test L ≈ [0.5890881713787511 0.7118839434795103]

    Q = eye_(2)
    R = eye_(1)
    L = dlqr(A,B,Q,R)
    @test L ≈ [0.5890881713787511 0.7118839434795103]

    B = [0;1]   # Note B is vector, B'B is scalar and INcompatible with matrix
    Q = eye_(2)
    R = eye_(1)
    @test_throws MethodError L ≈ dlqr(A,B,Q,R)
    #L ≈ [0.5890881713787511 0.7118839434795103]
end

end
