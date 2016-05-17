module TestMatrixComps
using CustomTest
using ControlSystems



A = [-0.21 0.2; 0.2 -0.21]
B = 0.01*eye(2)
C = eye(2)
D = 0
sys = ss(A,B,C,D)
sysr, G = balreal(sys)

@test_approx_eq gram(sysr, :c) G
@test_approx_eq gram(sysr, :o) G
@test_approx_eq sort(pole(sysr)) sort(pole(sys))

sysb,T = ControlSystems.balance_statespace(sys)
Ab,Bb,Cb,T = ControlSystems.balance_statespace(A,B,C)

@test_approx_eq sysb.A Ab

@test_approx_eq ControlSystems.balance_transform(A,B,C) ControlSystems.balance_transform(sys)


@test_approx_eq covar(sys, eye(2)) [0.002560975609756 0.002439024390244; 0.002439024390244 0.002560975609756]
D2 = eye(2)
sys2 = ss(A,B,C,D2, 1)
@test_approx_eq covar(sys2, eye(2)) [1.000110108378310 -0.000010098377310; -0.000010098377310 1.000110108378310]
#Direct term means infinite covariance
sys3 = ss(A,B,C,D2)
@test covar(sys3, eye(2)) == [Inf Inf; Inf Inf]

end
