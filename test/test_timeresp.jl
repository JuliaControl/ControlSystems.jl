import OrdinaryDiffEq

@testset "test_timeresp" begin

A = [0 1; 0 0]
B = [0;1]
C = [1 0]
sys = ss(A,B,C,0)
Q = I # eye(2)
R = I # eye(1)
L = lqr(sys,Q,R)

u(x,i) = -L*x # Form control law
t=0:0.1:50
x0 = [1.,0]
y, t, x, uout = lsim(sys,u,t,x0=x0) # Continuous time

th = 1e-6
@test sum(abs.(x[end,:])) < th

y, t, x, uout = lsim(c2d(sys,0.1)[1],u,t,x0=x0) # Discrete time
@test sum(abs.(x[end,:])) < th

#Do a manual simulation with uout
ym, tm, xm = lsim(sys, uout, t, x0=x0)
@test y ≈ ym
@test x ≈ xm

# Now compare to closed loop
# Discretization is needed before feedback
sysd = c2d(sys, 0.1)[1]
# Create the closed loop system
sysdfb = ss(sysd.A-sysd.B*L, sysd.B, sysd.C, sysd.D, 0.1)
#Simulate without input
yd, td, xd = lsim(sysdfb, zeros(501), t, x0=x0)
@test y ≈ yd
@test x ≈ xd

# Error for nonuniformly spaced vector
@test_throws ErrorException lsim(sys, [1, 2, 3, 4], [0, 1, 1, 2])

# Error for mismatch of sample rates of discrete-system and signal
@test_throws ErrorException lsim(sysd, [1, 2, 3, 4], 0:0.2:0.6)

# Test for problem with broadcast dimensions
@test lsim(sys, zeros(5), 0:0.2:0.8)[1][:] == zeros(5)


# lsim for Int system with Float64 input (regression test for #264)
@test lsim(ss(1,1,1,1,1), ones(5), 0:4)[1][:] == 1:5

# Various combinations of BigFloat
@test lsim(ss(big.(1),1,1,1,1), ones(5), 0:4)[1][:] == 1:5
@test lsim(ss(big.(1),1,1,1,1), big.(ones(5)), 0:4)[1][:] == 1:5
@test lsim(ss(1,1,1,1,1), big.(ones(5)), 0:4)[1][:] == 1:5

# Tests for HeteroStateSpace
@test lsim(HeteroStateSpace(big.(1.0),1,1,1,1), ones(5), 0:4)[1][:] == 1:5


# lsim for discrete-time complex-coefficient systems

# Complex system, real input signal
G1 = ss(1.0im, 1)
@test lsim(G1, [1, 2], 0:1)[1][:] == 1.0im*[1, 2]
@test impulse(G1, 4)[1][:] == [1.0im; zeros(4)]
@test step(G1, 4)[1][:] == fill(1.0im, 5)

G2 = ss(1.0im, 1, 1, 0, 1)
@test lsim(G2, ones(3), 0:2)[1][:] == [0.0, 1, 1 + im]
@test impulse(G2, 4)[1][:] == [0.0, 1, im, -1, -im]
@test step(G2, 4)[1][:] == [0.0, 1, 1+im, 1+im-1, 0]

# Real system, complex input signal
@test lsim(ss(1, 1), [1.0, 1.0im], 0:1)[1][:] == [1.0, 1.0im]
@test lsim(ss(1.0, 1), [1.0, 1.0im], 0:1)[1][:] == [1.0, 1.0im]
@test lsim(ss(1, 1, 1, 0, 1), [1.0, 1im, 1], 0:2)[1][:] == [0.0, 1, 1 + im]

# Complex system, complex input signal
@test lsim(ss(1.0im, 1), [1im, 2], 0:1)[1][:] == [-1, 2im]
@test lsim(ss(1.0im, 1, 1, 0, 1), [1.0, 1im, 1], 0:2)[1][:] == [0.0, 1, 2im]


# Test that the discrete lsim accepts u function that returns scalar
L = lqr(sysd,Q,R)
u(x,i) = -L*x
yd, td, xd = lsim(sysd, u, t, x0=x0)
@test norm(y - yd)/norm(y) < 0.05 # Since the cost matrices are not discretized, these will differ a bit
@test norm(x - xd)/norm(x) < 0.05

# Test lsim with default time vector
uv = randn(length(t))
y,t = lsim(c2d(sys,0.1)[1],uv,t,x0=x0)
yd,td = lsim(c2d(sys,0.1)[1],uv,x0=x0)
@test yd == y
@test td == t

#Test step and impulse Continuous
t0 = 0:0.05:2
systf = [tf(1,[1,1]) 0; 0 tf([1,-1],[1,2,1])]
sysss = ss([-1 0 0; 0 -2 -1; 0 1 0], [1 0; 0 1; 0 0], [1 0 0; 0 1 -1], 0)
yreal = zeros(length(t0), 2, 2)
xreal = zeros(length(t0), 3, 2)

#Step tf
y, t, x = step(systf, t0)
yreal[:,1,1] = 1 .- exp.(-t)
yreal[:,2,2] = -1 .+ exp.(-t)+2*exp.(-t).*t
@test y ≈ yreal atol=1e-4
#Step ss
y, t, x = step(sysss, t)
@test y ≈ yreal atol=1e-4
xreal[:,1,1] = yreal[:,1,1]
xreal[:,2,2] = exp.(-t).*t
xreal[:,3,2] = exp.(-t).*(-t .- 1) .+ 1
@test x ≈ xreal atol=1e-5

#Impulse tf
y, t, x = impulse(systf, t)
yreal[:,1,1] = exp.(-t)
yreal[:,2,2] = exp.(-t).*(1 .- 2*t)
@test y ≈ yreal atol=1e-2
#Impulse ss
y, t, x = impulse(sysss, t)
@test y ≈ yreal atol=1e-2
xreal[:,1,1] = yreal[:,1,1]
xreal[:,2,2] = -exp.(-t).*t + exp.(-t)
xreal[:,3,2] = exp.(-t).*t
@test x ≈ xreal atol=1e-1

#Test step and impulse Discrete
t0 = 0:0.05:2
systf = [tf(1,[1,1]) 0; 0 tf([1,-1],[1,2,1])]
sysss = ss([-1 0 0; 0 -2 -1; 0 1 0], [1 0; 0 1; 0 0], [1 0 0; 0 1 -1], 0)
yreal = zeros(length(t0), 2, 2)
xreal = zeros(length(t0), 3, 2)

#Step tf
y, t, x = step(systf, t0, method=:zoh)
yreal[:,1,1] = 1 .- exp.(-t)
yreal[:,2,2] = -1 .+ exp.(-t) + 2*exp.(-t).*t
@test y ≈ yreal atol=1e-14
#Step ss
y, t, x = step(sysss, t, method=:zoh)
@test y ≈ yreal atol=1e-13
xreal[:,1,1] = yreal[:,1,1]
xreal[:,2,2] = exp.(-t).*t
xreal[:,3,2] = exp.(-t).*(-t .- 1) .+ 1
@test x ≈ xreal atol=1e-14

#Impulse tf
y, t, x = impulse(systf, t, method=:zoh)
yreal[:,1,1] = exp.(-t)
yreal[:,2,2] = exp.(-t).*(1 .- 2*t)
@test y ≈ yreal atol=1e-14
#Impulse ss
y, t, x = impulse(1.0sysss, t, method=:zoh)
@test y ≈ yreal atol=1e-13
xreal[:,1,1] = yreal[:,1,1]
xreal[:,2,2] = -exp.(-t).*t + exp.(-t)
xreal[:,3,2] = exp.(-t).*t
@test x ≈ xreal atol=1e-13


#Step response of discrete system with specified final time
G = tf([1], [1; zeros(3)], 1)
y, t2, x = step(G, 10)
@test y ≈ [zeros(3); ones(8)] atol=1e-5
@test t2 == 0:1:10 # isapprox is broken for ranges (julia 1.3.1)

#Impulse response of discrete system to final time that is not mulitple of the sample time
G = tf([1], [1; zeros(3)], 0.3)
y, t2, x = step(G, 2)
@test y ≈ [zeros(3); ones(4)] atol=1e-5
@test t2 ≈ 0:0.3:1.8 atol=1e-5

#Make sure t was never changed
@test t0 == t








@testset "Simulators" begin


h              = 0.1
Tf             = 20
t              = 0:h:Tf
P              = ss(tf(1,[2,1])^2)
reference(x,t) = [1.]
s              = Simulator(P, reference)
x0             = [0.,0]
tspan          = (0.0,Tf)
sol            = solve(s, x0, tspan, OrdinaryDiffEq.Tsit5())
@test step(P,t)[3] ≈ reduce(hcat,sol.(t))'
end

end
