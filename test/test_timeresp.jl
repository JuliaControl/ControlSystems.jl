import OrdinaryDiffEq
using OrdinaryDiffEq: solve
using ControlSystems
using Plots

A = [0 1; 0 0]
B = [0;1]
C = [1 0]
sys = ss(A,B,C,0)
Q = I # eye(2)
R = I # eye(1)
L = lqr(sys,Q,R)

u1(x,i) = -L*x # Form control law
t=0:0.1:50
x0 = [1.,0]
res = lsim(sys,u1,t,x0=x0) # Continuous time
y, t, x, uout = res

th = 1e-6
@test sum(abs.(x[:,end])) < th


# Error for wrong size of u
@test_throws ErrorException lsim(sys, (x,t)->[sin(t), cos(t)], 2; x0=[0, 1])

# Restrict u to being vectors
@test_throws ErrorException lsim(sys, (x,t)->[sin(t)]', 0:0.1:1)
@test_throws ErrorException lsim(sys, (x,t)->[sin(t)]', 2)
@test_throws ErrorException lsim(sys, (x,t)->sin(t), 2)

# Error for wrong size of x0
@test_throws ErrorException lsim(sys, (x,t)->[sin(t)], 2; x0=[0, 1, 1])

# Error for nonuniformly spaced vector
@test_throws ErrorException lsim(sys, [1 2 3 4], [0, 1, 1, 2])

# Test for problem with broadcast dimensions
@test lsim(sys, zeros(1, 5), 0:0.2:0.8)[1][:] == zeros(5)

#Test step and impulse Continuous
t0 = 0:0.05:2
systf = [tf(1,[1,1]) 0; 0 tf([1,-1],[1,2,1])]
sysss = ss([-1 0 0; 0 -2 -1; 0 1 0], [1 0; 0 1; 0 0], [1 0 0; 0 1 -1], 0)
yreal = zeros(2, length(t0), 2)
xreal = zeros(3, length(t0), 2)

#Step tf
y, t, x = step(systf, t0)
yreal[1,:,1] = 1 .- exp.(-t)
yreal[2,:,2] = -1 .+ exp.(-t)+2*exp.(-t).*t
@test y ≈ yreal atol=1e-4
#Step ss
y, t, x = step(sysss, t)
@test y ≈ yreal atol=1e-4
xreal[1,:,1] = yreal[1,:,1]
xreal[2,:,2] = exp.(-t).*t
xreal[3,:,2] = exp.(-t).*(-t .- 1) .+ 1
@test x ≈ xreal atol=1e-5

#Impulse tf
y, t, x = impulse(systf, t)
yreal[1,:,1] = exp.(-t)
yreal[2,:,2] = exp.(-t).*(1 .- 2*t)
@test y ≈ yreal atol=1e-2
#Impulse ss
y, t, x = impulse(sysss, t)
@test y ≈ yreal atol=1e-2
xreal[1,:,1] = yreal[1,:,1]
xreal[2,:,2] = -exp.(-t).*t + exp.(-t)
xreal[3,:,2] = exp.(-t).*t
@test x ≈ xreal atol=1e-1

#Step tf
y, t, x = step(systf, t0, method=:zoh)
yreal[1,:,1] = 1 .- exp.(-t)
yreal[2,:,2] = -1 .+ exp.(-t) + 2*exp.(-t).*t
@test y ≈ yreal atol=1e-14
#Step ss
y, t, x = step(sysss, t, method=:zoh)
@test y ≈ yreal atol=1e-13
xreal[1,:,1] = yreal[1,:,1]
xreal[2,:,2] = exp.(-t).*t
xreal[3,:,2] = exp.(-t).*(-t .- 1) .+ 1
@test x ≈ xreal atol=1e-14

#Impulse tf
y, t, x = impulse(systf, t, method=:zoh)
yreal[1,:,1] = exp.(-t)
yreal[2,:,2] = exp.(-t).*(1 .- 2*t)
@test y ≈ yreal atol=1e-14
#Impulse ss
y, t, x = impulse(1.0sysss, t, method=:zoh)
@test y ≈ yreal atol=1e-13
xreal[1,:,1] = yreal[1,:,1]
xreal[2,:,2] = -exp.(-t).*t + exp.(-t)
xreal[3,:,2] = exp.(-t).*t
@test x ≈ xreal atol=1e-13


@testset "Simulators" begin
    Ts             = 0.1
    tfinal         = 20
    t              = 0:Ts:tfinal
    P              = ss(tf(1,[2,1])^2)
    reference(x,t) = [1.]
    s              = Simulator(P, reference)
    x0             = [0.,0]
    tspan          = (0.0,tfinal)
    sol            = solve(s, x0, tspan, OrdinaryDiffEq.Tsit5())
    @test step(P,t)[3] ≈ reduce(hcat,sol.(t)) rtol=1e-4
end

# test fwd euler
C_111 = ss([-5], [2], [3], [0])
Cd = c2d(C_111, 0.001, :fwdeuler)
t = 0:0.001:2
y,_ = step(C_111, t)
yd,_ = step(Cd, t)
@test norm(yd-y) / norm(y) < 1e-3


## Test time scaling
for balanced in [true, false]
    sys = ssrand(1,1,5);
    t = 0:0.1:50
    a = 10
    res1 = step(sys, t)
    sys2 = time_scale(sys, a; balanced)
    res2 = step(sys2, t ./ a)
    @test res1.y ≈ res2.y rtol=1e-2 atol=1e-2
end


## stepinfo
G = tf(1, [1, 1, 1])
res = step(G, 20)
si = stepinfo(res)
show(si)
@test si.settlingtime ≈ 8.134 rtol=0.01
@test si.peak ≈ 1.163 rtol=0.01
@test si.peaktime ≈ 3.652 rtol=0.01
@test si.overshoot ≈ 16.295 rtol=0.01
@test si.undershoot == 0
@test si.risetime ≈ 1.660 rtol=0.01
@test si.yf ≈ 1.0 rtol=0.01
@test_nowarn plot(si)

# test reverse direction
G = tf(-1, [1, 1, 1])
res = step(G, 20)
si2 = stepinfo(res)
@test si2.settlingtime ≈ si.settlingtime rtol=0.01
@test si2.peak ≈ -si.peak rtol=0.01
@test si2.peaktime ≈ si.peaktime rtol=0.01
@test si2.overshoot ≈ si.overshoot rtol=0.01
@test si.undershoot == 0
@test si2.risetime ≈ si.risetime rtol=0.01
@test si2.yf ≈ -si.yf rtol=0.01
@test_nowarn plot(si2)

# very poorly damped
G = tf(1, [1, 0.01, 1])
res = step(G, 200)
@test_logs (:warn, "System might not have settled within the simulation time") stepinfo(res, yf=1)
si = stepinfo(res, yf=1) # Provide final value manually
@test_nowarn plot(si)

# non minimum phase
G = tf([-1, 1], [1, 1, 1])
res = step(G, 20)
si = stepinfo(res)
test si.undershoot ≈ 27.98 rtol=0.01
plot(si)