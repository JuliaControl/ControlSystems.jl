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

u(t,x) = -L*x # Form control law
t=0:0.1:50
x0 = [1,0]
y, t, x, uout = lsim(sys,u,t,x0)
@test sum(abs(x[end,:])) < eps()
# plot(t,x, lab=["Position", "Velocity"]', xlabel="Time [s]")

end
