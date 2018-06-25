@testset "test_matrix_comps" begin
A = [-0.21 0.2; 0.2 -0.21]
B = 0.01*eye(2)
C = eye(2)
D = 0
sys = ss(A,B,C,D)
sysr, G = balreal(sys)

@test gram(sysr, :c) ≈ G
@test gram(sysr, :o) ≈ G
@test sort(pole(sysr)) ≈ sort(pole(sys))

sysb,T = ControlSystems.balance_statespace(sys)
Ab,Bb,Cb,T = ControlSystems.balance_statespace(A,B,C)

@test Ab*T ≈ T*A
@test Bb ≈ T*B
@test Cb*T ≈ C

@test sysb.A ≈ Ab

@test ControlSystems.balance_transform(A,B,C) ≈ ControlSystems.balance_transform(sys)

W = eye(2)
@test covar(sys, W) ≈ [0.002560975609756 0.002439024390244; 0.002439024390244 0.002560975609756]
D2 = eye(2)
@test covar(ss(A,B,C,D2, 1), W) ≈ [1.000110108378310 -0.000010098377310; -0.000010098377310 1.000110108378310]
# Direct term means infinite covariance
@test covar(ss(A,B,C,D2), W) ≈ [Inf Inf; Inf Inf]

# No noise on second output should give finite variance
@test covar(ss(A,B,C,[1 0; 0 0]), W) ≈ [Inf Inf; Inf 0.002560975609756]

# Unstable system has inf covar
@test covar(ss(eye(2),B,C,0), W) == [Inf Inf; Inf Inf]

# Discrete system can have direct term
@test covar(ss(A,B,C,D2,0.1),W) ≈ [1.00011010837831 -1.0098377309782909e-5; -1.0098377309782909e-5 1.00011010837831]

# TODO test in Julia 0.7 to see if supported
# # Test special matrices
# As = sparse(A)
# Bs = sparse(B)
# Cs = sparse(C)
# Asb,Bsb,Csb,Ts = ControlSystems.balance_statespace(As,Bs,Cs) #Error no LAPACK function
#
# @test Abs*Ts ≈ Ts*As
# @test Bbs ≈ Ts*Bs
# @test Cbs*Ts ≈ Cs

# Test special values
Ar = rationalize.(A)
Br = rationalize.(B)
Cr = rationalize.(C)
Arb,Brb,Crb,Tr = ControlSystems.balance_statespace(Ar,Br,Cr)

@test Arb*Tr ≈ Tr*Ar
@test Brb ≈ Tr*Br
@test Crb*Tr ≈ Cr

end
