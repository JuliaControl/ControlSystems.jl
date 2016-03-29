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

end
