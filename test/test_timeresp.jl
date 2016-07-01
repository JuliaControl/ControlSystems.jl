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


#Test step and impulse
t0 = 0:0.05:2
systf = [tf(1,[1,1]) 0; 0 tf([1,-1],[1,2,1])]
sysss = ss([-1 0 0; 0 -2 -1; 0 1 0], [1 0; 0 1; 0 0], [1 0 0; 0 1 -1], 0)
yreal = zeros(length(t0), 2, 2)
xreal = zeros(length(t0), 3, 2)

#Step tf
y, t, x = step(systf, t0)
yreal[:,1,1] = 1-exp(-t)
yreal[:,2,2] = -1+exp(-t)+2*exp(-t).*t
@test_approx_eq_eps y yreal 1e-14
#Step ss
y, t, x = step(sysss, t)
@test_approx_eq_eps y yreal 1e-14
xreal[:,1,1] = yreal[:,1,1]
xreal[:,2,2] = exp(-t).*t
xreal[:,3,2] = exp(-t).*(-t-1) + 1
@test_approx_eq_eps x xreal 1e-14

#Impulse tf
y, t, x = impulse(systf, t)
yreal[:,1,1] = exp(-t)
yreal[:,2,2] = exp(-t).*(1 - 2.*t)
@test_approx_eq_eps y yreal 1e-14
#Impulse ss
y, t, x = impulse(sysss, t)
@test_approx_eq_eps y yreal 1e-14
xreal[:,1,1] = yreal[:,1,1]
xreal[:,2,2] = -exp(-t).*t + exp(-t)
xreal[:,3,2] = exp(-t).*t
@test_approx_eq_eps x xreal 1e-14


#Step response of discrete system with specified final time
G = tf([1], [1; zeros(3)], 1)
y, t2, x = step(G, 10)
@test_approx_eq_eps y [zeros(3); ones(8)] 1e-14
@test_approx_eq_eps t2 0:1:10 1e-14

#Impulse response of discrete system to final time that is not mulitple of the sample time 
G = tf([1], [1; zeros(3)], 0.3)
y, t2, x = step(G, 2)
@test_approx_eq_eps y [zeros(3); ones(4)] 1e-14
@test_approx_eq_eps t2 0:0.3:1.8 1e-14



#Make sure t was never changed
@test t0 == t
end
