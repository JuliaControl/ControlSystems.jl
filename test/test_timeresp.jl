module TestTimeResp
using CustomTest
using ControlSystems

A = [0 1; 0 0]
B = [0;1]
C = [1 0]
sys = ss(A,B,C,0)
Q = eye(2)
R = eye(1)
L = lqr(sys,Q,R)

u(i,x) = -L*x # Form control law
t=0:0.1:50
x0 = [1,0]
y, t, x, uout = lsim(sys,u,t,x0)
@test sum(abs(x[end,:])) < eps()

#Do a manual simulation with uout
ym, tm, xm = lsim(sys, uout, t, x0)
@test_approx_eq y ym
@test_approx_eq x xm

# Now compare to closed loop
# Discretization is needed before feedback
sysd = c2d(sys, 0.1)[1]
# Create the closed loop system
sysdfb = ss(sysd.A-sysd.B*L, sysd.B, sysd.C, sysd.D, 0.1)
#Simulate without input
yd, td, xd = lsim(sysdfb, zeros(501), t, x0)

@test_approx_eq y yd
@test_approx_eq x xd

end
