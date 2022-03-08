@testset "test_matrix_comps" begin
A = [-0.21 0.2; 0.2 -0.21]
B = 0.01*[1 0; 0 1]
C = [1 0; 0 1]
D = 0
sys = ss(A,B,C,D)
sysr, G = balreal(sys)

@test gram(sysr, :c) ≈ G
@test gram(sysr, :o) ≈ G
@test sort(poles(sysr)) ≈ sort(poles(sys))

sysb,T = ControlSystems.balance_statespace(sys)
@test similarity_transform(sysb, T) ≈ sys
Ab,Bb,Cb,T = ControlSystems.balance_statespace(A,B,C)

@test Ab*T ≈ T*A
@test Bb ≈ T*B
@test Cb*T ≈ C

@test sysb.A ≈ Ab
@test similarity_transform(sysb, T) ≈ sys


U = svd(randn(2,2)).U
syst = similarity_transform(sys, U, unitary = true)
Ab,Bb,Cb,Db = ssdata(syst)
@test Ab ≈ U'A*U
@test Bb ≈ U'B
@test Cb ≈ C*U


@test ControlSystems.balance_transform(A,B,C) ≈ ControlSystems.balance_transform(sys)

W = [1 0; 0 1]
@test covar(sys, W) ≈ [0.002560975609756 0.002439024390244; 0.002439024390244 0.002560975609756]
D2 = [1 0; 0 1]
@test covar(ss(A,B,C,D2, 1), W) ≈ [1.000110108378310 -0.000010098377310; -0.000010098377310 1.000110108378310]
# Direct term means infinite covariance
@test covar(ss(A,B,C,D2), W) ≈ [Inf Inf; Inf Inf]

# No noise on second output should give finite variance
@test covar(ss(A,B,C,[1 0; 0 0]), W) ≈ [Inf Inf; Inf 0.002560975609756]

# Unstable system has inf covar
@test covar(ss([1 0; 0 1],B,C,0), W) == [Inf Inf; Inf Inf]

# Discrete system can have direct term
@test covar(ss(A,B,C,D2,0.1),W) ≈ [1.00011010837831 -1.0098377309782909e-5; -1.0098377309782909e-5 1.00011010837831]

# Static systems
R = fill(1.0,1,1)
@test covar(ss(2.0), R) == fill(Inf, 1, 1)
@test covar(ss(2.0, 0.2), R) == fill(4.0, 1, 1)


# TODO test in Julia 0.7 to see if supported
# # Test special matrices
As = sparse(A)
Bs = sparse(B)
Cs = sparse(C)
@test_logs (:warn, "Unable to balance state-space, returning original system") ControlSystems.balance_statespace(As,Bs,Cs)
#
# @test Abs*Ts ≈ Ts*As
# @test Bbs ≈ Ts*Bs
# @test Cbs*Ts ≈ Cs

# Test special values
Ar = rationalize.(A)
Br = rationalize.(B)
Cr = rationalize.(Float64.(C))    # When did ever rationalize work on Int?
Arb,Brb,Crb,Tr = ControlSystems.balance_statespace(Ar,Br,Cr)

@test Arb*Tr ≈ Tr*Ar
@test Brb ≈ Tr*Br
@test Crb*Tr ≈ Cr

Tr = randn(2,2)
syst = similarity_transform(sys, Tr)
@test sys.A*Tr ≈ Tr*syst.A
@test sys.B ≈ Tr*syst.B
@test sys.C*Tr ≈ syst.C


sys = ss([1 0.1; 0 1], ones(2), [1. 0], 0)
sysi = ControlSystems.innovation_form(sys, I, I)
@test sysi.A ≈ sysi.A
@test sysi.B ≈ [4.415675759647131
 48.334204475215365]

 sysi = innovation_form(sys)
 @test sysi.B ≈ [4.415675759647131
 48.334204475215365]

 sysi = innovation_form(sys, R2 = 2I)
 @test sysi.B ≈ [4.225661436353894
 44.52445850991302]

 sysi = innovation_form(sys, R1 = 2I)
 @test sysi.B ≈ [4.734159731874057
 54.719744515739514]

 sysi = ControlSystems.innovation_form(sys, R1=2I, R2=2I)
 @test sysi.B ≈ [4.415675759647131
 48.334204475215365]

# Test with noise filters
sysw = ss([0.5 0.1; 0 0.5], [0,1], eye_(2), 0, 1)
sysi = ControlSystems.innovation_form(sys, sysw=sysw)
@test sysi.A ≈ sys.A
@test sysi.B ≈ [4.01361818808572
 40.26132476965486]


# Test observer_predictor
sysp = ControlSystems.observer_predictor(sys, I(2), I(1))
K = kalman(sys, I(2), I(1))
@test sysp.A == sys.A-K*sys.C
@test sysp.B == [sys.B-K*sys.D K]

# test longer prediction horizons

# With K=0, y should have no influence
sys = ssrand(1,1,2, Ts=1)
u = randn(1,5)
K = zeros(2,1)
sysp = observer_predictor(sys, K; h=2)

# test with two different random outputs
@test lsim(sysp, [u; randn(1,5)]).y == lsim(sysp, [u; randn(1,5)]).y

# test that it's the same as simulating the system by itself
@test lsim(sysp, [u; randn(1,5)]).y == lsim(sys, u).y


# With K != 0 but u = 0
K = randn(2,1)
h = 3
sysp = observer_predictor(sys, K; h)
u = zeros(1, 5)
y = zeros(1, 5)
y[1] = 1

yh = lsim(sysp, [u; y]).y

# The first h outputs are all zero
@test maximum(abs, yh[1:h]) < 10eps()
@test abs(yh[h+1]) > eps()

##
A = [
    -0.6437   -0.5055   -0.3211  -0.03438
    -0.165    0.7435   -0.4341   -0.2137
    -0.07843    0.1487    0.4797   -0.7199
    0.5509   -0.3798    -0.175   -0.4614]

B = [
    -0.001511   -0.03363
    0.01829   -0.01114
    -0.03234   -0.01333
    0.01556    -0.0131]

C = [-33.01   32.88   23.15  -5.471]

D = [0   0]

K = [
    -0.005343
    0.003431
    0.0139
    0.008578;;]

sys = ss(A,B,C,D,1)
sysp = observer_predictor(sys, K; h=3)
u = [0.7378277077784556 -0.9823648579212658 0.16545080647905613 -0.04218269410737019 2.0392261068878264 1.1557293975007483 -0.38811443803683415 2.0733577162855688 -0.03625300688766673 -0.12317844047992449; 0.14799957126508187 0.4859599561753988 -1.07200839703822 -0.8326580050678177 0.24617707129291685 0.6912641690407068 0.3761998517214711 -1.604656172130068 -0.747064482343228 -0.974687925632395]
y = [-0.842859151936594 -0.19786693123277427 0.6436197988164484 0.053887981499727844 1.141282640261743 -0.9797169106947525 -0.19891450532826468 1.2050598007441486 0.9346495338941443 0.48283808371749715]

res = lsim(sysp, [u; y])
@test res.y ≈ [0.0 -0.05966333850248531 0.503489609536266 -0.20592242042526876 0.5894787365287978 -0.6386410889293671 0.45779714753135586 2.278268660780028 2.4541372811479922 4.849094227405488]



# Test observer_controller
sys = ssrand(2,3,4)
Q1 = I(4)
Q2 = I(3)
R1 = I(4)
R2 = I(2)
L = lqr(sys, Q1, Q2)
K = kalman(sys, R1, R2)
cont = observer_controller(sys, L, K)
syscl = feedback(sys, cont)

pcl = poles(syscl)
A,B,C,D = ssdata(sys)
allpoles = [
    eigvals(A-B*L)
    eigvals(A-K*C)
]
@test sort(pcl, by=LinearAlgebra.eigsortby) ≈ sort(allpoles, by=LinearAlgebra.eigsortby) 
@test cont.B == K

end 
