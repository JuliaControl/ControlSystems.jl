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

@testset "placemimo" begin
Random.seed!(0)
A = randn(3, 3)
B = randn(3, 4)

p = [-3.0, -2, 1.3]
for m in 1:4 # Test for m=1, 1<m<n, m=n, m>n
    K = ControlSystems.placemimo(A, B[:, 1:m], p)
    @test eigvals(A-B[:, 1:m]*K) ≈ p
end

p = [-3.0, -2, -2]
K = ControlSystems.placemimo(A, B[:, 1:2], p)
@test eigvals(A-B[:, 1:2]*K) ≈ p
K = ControlSystems.placemimo(A, B[:, 1:3], p)
@test eigvals(A-B[:, 1:3]*K) ≈ p

p = [-2.0, -1-im, -1+im]
K = ControlSystems.placemimo(A, B[:, 1:1], p)
@test eigvals(A-B[:, 1:1]*K) ≈ p
K = ControlSystems.placemimo(A, B[:, 1:2], p)
@test eigvals(A-B[:, 1:2]*K) ≈ p

# Should not be able to place if there is a pole with multiplicity higher than m given that (A, B) is controllable
p = [-3.0, -2, -2]
@test_throws AssertionError ControlSystems.placemimo(A, B[:, 1:1], p)

# Not sure how interesting this is
p = [-3.0, -2, 1.3]
@test_throws ErrorException ControlSystems.placemimo(A, B[:, 1:2], p, max_iter=1)

p = [-3.0, -2, -2]
@inferred ControlSystems.placemimo(A, B[:, 1:2], p)
p = [-2.0, -1-im, -1+im]
@inferred ControlSystems.placemimo(A, B[:, 1:2], p)
end
# test where B is not full rank but it is still possible? should maybe not work according to article? test throws?
# test where A, B not controllable? but controllable mode is in poles, should work
# test where B = nxm and m < multiplicity(pole_i) which should not work
# testa vilken som är effektivast för single column B
# en test throws @test_throws som visar bra felmeddelande

@testset "LQR" begin
    h = 0.1
    A = [1 h; 0 1]
    B = [0;1] # Note B is vector, B'B is scalar, but compatible with I
    C = [1 0]
    Q = I
    R = I
    L = dlqr(A,B,Q,R)
    @test L ≈ [0.5890881713787511 0.7118839434795103]
    sys = ss(A,B,C,0,h)
    L = lqr(sys, Q, R)
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
