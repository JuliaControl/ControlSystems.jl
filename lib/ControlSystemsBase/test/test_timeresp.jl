using ControlSystemsBase

@testset "test_timeresp" begin

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
@test_throws MethodError lsim(sys,u1,t,x0=x0) # Continuous time simulation not loaded


th = 1e-6

sysd = c2d(sys,0.1)
res = lsim(sysd,u1,t,x0=x0) # Discrete time
y, t, x, uout = res
@test sum(abs.(x[:,end])) < th

## Use workspace method directly
lsws = LsimWorkspace(sysd, length(t))
res2 = lsim!(lsws, sysd, u1, t; x0 = copy(x0))
@test res2.sys === res.sys
@test res2.y == res.y
@test res2.x == res.x
@test res2.u == res.u

#Do a manual simulation with uout
resm = lsim(sys, uout, t, x0=x0)
ym, tm, xm = resm
@test y ≈ ym
@test x ≈ xm

# Discrete-time
resm = lsim(sysd, uout, t, x0=x0)
ym, tm, xm = resm
@test y ≈ ym
@test x ≈ xm

#Do a manual simulation with uout and workspace
res2 = lsim!(lsws, sysd, uout, t; x0 = copy(x0))
@test res2.sys === resm.sys
@test res2.y == resm.y
@test res2.x == resm.x
@test res2.u == resm.u

@test lsim!(lsws, sysd, uout, t; x0 = copy(x0)).t == res2.t

# Now compare to closed loop
# Discretization is needed before feedback
# Create the closed loop system
sysdfb = ss(sysd.A-sysd.B*L, sysd.B, sysd.C, sysd.D, 0.1)
#Simulate without input
yd, td, xd = lsim(sysdfb, zeros(1, 501), t, x0=x0)
@test y ≈ yd
@test x ≈ xd

# Error for nonuniformly spaced vector
@test_throws ErrorException lsim(sys, [1 2 3 4], [0, 1, 1, 2])

# Error for mismatch of sample rates of discrete-system and signal
@test_throws ErrorException lsim(sysd, [1 2 3 4], 0:0.2:0.6)

# Warning for wrong input dimension
@test_logs (:warn, r"should be a matrix") lsim(sys, [1; 2; 3; 4;;], 1:4)
@test_throws TypeError lsim(sys, [1; 2; 3; 4], 1:4)
@test_logs (:warn, r"row-vector") lsim(sysd, [1, 2, 3, 4])

# Test for problem with broadcast dimensions
@test lsim(sys, zeros(1, 5), 0:0.2:0.8)[1][:] == zeros(5)


# lsim for Int system with Float64 input (regression test for #264)
@test lsim(ss(1,1,1,1,1), ones(1, 5), 0:4)[1][:] == 1:5

# Various combinations of BigFloat
@test lsim(ss(big.(1),1,1,1,1), ones(1, 5), 0:4)[1][:] == 1:5
@test lsim(ss(big.(1),1,1,1,1), big.(ones(1, 5)), 0:4)[1][:] == 1:5
@test lsim(ss(1,1,1,1,1), big.(ones(1, 5)), 0:4)[1][:] == 1:5

# Tests for HeteroStateSpace
sysb = HeteroStateSpace(big.(1.0),1,1,1,1)
u = ones(1, 5)
@test lsim(sysb, u, 0:4)[1][:] == 1:5
ws = LsimWorkspace(sysb, u)
@test lsim!(ws, sysb, u, 0:4)[1][:] == 1:5

# test explicit size and type constructor
ws2 = LsimWorkspace{Float32}(1,2,3,10)
@test size(ws2.y) == (1,10)
@test size(ws2.u) == (2,10)
@test size(ws2.x) == (3,10)
@test ws2.x isa Array{Float32, 2}

# lsim for discrete-time complex-coefficient systems

# Complex system, real input signal
G1 = ss(1.0im, 1)
@test lsim(G1, [1 2], 0:1)[1][:] == 1.0im*[1, 2]
@test impulse(G1, 4)[1][:] == [1.0im; zeros(4)]
@test step(G1, 4)[1][:] == fill(1.0im, 5)

G2 = ss(1.0im, 1, 1, 0, 1)
@test lsim(G2, ones(1, 3), 0:2)[1][:] == [0.0, 1, 1 + im]
@test impulse(G2, 4)[1][:] == [0.0, 1, im, -1, -im]
@test step(G2, 4)[1][:] == [0.0, 1, 1+im, 1+im-1, 0]

# Real system, complex input signal
@test lsim(ss(1, 1), [1.0 1.0im], 0:1)[1][:] == [1.0, 1.0im]
@test lsim(ss(1.0, 1), [1.0 1.0im], 0:1)[1][:] == [1.0, 1.0im]
@test lsim(ss(1, 1, 1, 0, 1), [1.0 1im 1], 0:2)[1][:] == [0.0, 1, 1 + im]

# Complex system, complex input signal
@test lsim(ss(1.0im, 1), [1im 2], 0:1)[1][:] == [-1, 2im]
@test lsim(ss(1.0im, 1, 1, 0, 1), [1.0 1im 1], 0:2)[1][:] == [0.0, 1, 2im]


# Test that the discrete lsim accepts u function that returns scalar
L = lqr(sysd,Q,R)
u2(x,i) = -L*x
yd, td, xd = lsim(sysd, u2, t, x0=x0)
@test norm(y - yd)/norm(y) < 0.05 # Since the cost matrices are not discretized, these will differ a bit
@test norm(x - xd)/norm(x) < 0.05

# Test lsim with default time vector
uv = randn(1, length(t))
y,t = lsim(c2d(sys,0.1),uv,t,x0=x0)
yd,td = lsim(c2d(sys,0.1),uv,x0=x0)
@test yd == y
@test td == t


#Test step and impulse Discrete
t0 = 0:0.05:2
systf = [tf(1,[1,1]) 0; 0 tf([1,-1],[1,2,1])]
sysss = ss([-1 0 0; 0 -2 -1; 0 1 0], [1 0; 0 1; 0 0], [1 0 0; 0 1 -1], 0)
yreal = zeros(2, length(t0), 2)
xreal = zeros(3, length(t0), 2)

#Step tf
y, t, x = step(systf, t0, method=:zoh) # Method is :zoh so discretization is applied
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


#Step response of discrete system with specified final time
G = tf([1], [1; zeros(3)], 1)
y, t2, x = step(G, 10)
@test y ≈ [zeros(1, 3) ones(1, 8)] atol=1e-5
@test t2 == 0:1:10 # isapprox is broken for ranges (julia 1.3.1)

#Impulse response of discrete system to final time that is not mulitple of the sample time
G = tf([1], [1; zeros(3)], 0.3)
y, t2, x = step(G, 2)
@test y ≈ [zeros(1, 3) ones(1, 4)] atol=1e-5
@test t2 ≈ 0:0.3:1.8 atol=1e-5

#Make sure t was never changed
@test t0 == t


# Test _default_dt 
sysstab = tf(1, [1, 0.5])
sysstatic = zpk([], [], 1.0)
sysunstab = tf(1, [1, -1])
sysint = tf(1, [1, 0])
sysd = tf(1, [1, 1], 0.01)

sys2tau = tf(1, [0.1, 1]) * tf(1, [10, 1])

@test ControlSystemsBase._default_dt(sysstab) == 0.17
@test ControlSystemsBase._default_dt(sysstatic) == 0.05
@test ControlSystemsBase._default_dt(sysunstab) == 0.083
@test ControlSystemsBase._default_dt(sysint) == 0.05
@test ControlSystemsBase._default_dt(sysd) == 0.01
@test 10 < ControlSystemsBase._default_time_vector(sysstab)[end] < 100 # These are set loose to allow some future tweaking. 
@test 5 < ControlSystemsBase._default_time_vector(sysint)[end] < 15
@test 20 < ControlSystemsBase._default_time_vector(sys2tau)[end] < 100


# Test error hints
if VERSION >= v"1.7"
    # If ControlSystems is not loaded, the following should error
    G = ssrand(1,1,1)
    @test_throws "install and load ControlSystems.jl" lsim(G, (u,t)->[1], 10)
    @test_throws "install and load ControlSystems.jl" lsim(G*delay(1), (u,t)->[1], 10)
    @test_throws "step with continuous-time systems" step(G*delay(1), 10)
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
@test si.undershoot ≈ 27.98 rtol=0.01
plot(si)

end
