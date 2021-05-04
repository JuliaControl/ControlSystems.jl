import DelayDiffEq: MethodOfSteps, Tsit5

@testset "test_delay_system" begin

# For simplicity, equality of DelayLtiSystems are tested over a finite set of frequencies
ω = 0.0:8


@test typeof(promote(delay(0.2), ss(1))[1]) == DelayLtiSystem{Float64,Float64}

@test typeof(promote(delay(0.2), ss(1.0 + im))[1]) == DelayLtiSystem{Complex{Float64}, Float64}

if VERSION >= v"1.6.0-DEV.0"
    @test sprint(show, ss(1,1,1,1)*delay(1.0)) == "DelayLtiSystem{Float64, Float64}\n\nP: StateSpace{Continuous, Float64}\nA = \n 1.0\nB = \n 0.0  1.0\nC = \n 1.0\n 0.0\nD = \n 0.0  1.0\n 1.0  0.0\n\nContinuous-time state-space model\n\nDelays: [1.0]\n"
else
    @test sprint(show, ss(1,1,1,1)*delay(1.0)) == "DelayLtiSystem{Float64,Float64}\n\nP: StateSpace{Continuous,Float64}\nA = \n 1.0\nB = \n 0.0  1.0\nC = \n 1.0\n 0.0\nD = \n 0.0  1.0\n 1.0  0.0\n\nContinuous-time state-space model\n\nDelays: [1.0]\n"
end

# Extremely baseic tests
@test freqresp(delay(1), ω) ≈ reshape(exp.(-im*ω), length(ω), 1, 1) rtol=1e-15
@test freqresp(delay(2.5), ω)[:] ≈ exp.(-2.5im*ω) rtol=1e-15
@test freqresp(3.5*delay(2.5), ω)[:] ≈ 3.5*exp.(-2.5im*ω) rtol=1e-15
@test freqresp(delay(2.5)*1.5, ω)[:] ≈ exp.(-2.5im*ω)*1.5 rtol=1e-15

# Stritcly proper system
P1 = DelayLtiSystem(ss(-1.0, 1, 1, 0))
P1_fr = 1 ./ (im*ω .+ 1)
@test freqresp(P1, ω)[:] == P1_fr

# Not stritcly proper system
P2 = DelayLtiSystem(ss(-2.0, -1, 1, 1)) # (s+1)/(s+2)
P2_fr = (im*ω .+ 1) ./ (im*ω .+ 2)
@test freqresp(P2, ω)[:] ≈ P2_fr rtol=1e-15


## Addition
@test freqresp(1 + delay(1), ω)[:] ≈ 1 .+ exp.(-im*ω)
@test freqresp(P1 + delay(1), ω)[:] ≈ P1_fr .+ exp.(-im*ω)

# Substraction
@test freqresp(P1 - delay(1), ω)[:] ≈ P1_fr .- exp.(-im*ω)

## Multiplication by scalar
@test freqresp(2*delay(1), ω)[:] ≈ 2*exp.(-im*ω)
@test freqresp(2*(P1*delay(1)), ω)[:] ≈ 2*P1_fr.*exp.(-im*ω)

## Multiplication
@test freqresp(P1 * delay(1), ω)[:] ≈ P1_fr .* exp.(-im*ω) rtol=1e-15
@test freqresp(delay(1) * P1, ω)[:] ≈ P1_fr .* exp.(-im*ω) rtol=1e-15

@test freqresp(P2 * delay(1), ω)[:] ≈ P2_fr .* exp.(-im*ω) rtol=1e-15
@test freqresp(delay(1) * P2, ω)[:] ≈ P2_fr .* exp.(-im*ω) rtol=1e-15


# Equality
@test P1 == deepcopy(P1)
@test P1 != deepcopy(P2)

# evalfr
s_vec = [0, 1im, 1, 1 + 1im]
@test [evalfr(delay(2), s)[1] for s in s_vec] ≈ [exp(-2*s) for s in s_vec] rtol=1e-16

## Feedback
# The first tests don't include delays, but the linear system is of DelayLtiForm type
# (very simple system so easy to troubleshoot)
@test freqresp(feedback(1.0, P1), ω)[:] ≈ 1 ./ (1 .+  P1_fr) rtol=1e-15
@test freqresp(feedback(P1, 1.0), ω)[:] ≈ P1_fr ./ (1 .+  P1_fr) rtol=1e-15
@test freqresp(feedback(P1, P1), ω)[:] ≈ P1_fr ./ (1 .+ P1_fr .* P1_fr) rtol=1e-15

@test freqresp(feedback(1.0, DelayLtiSystem(ss(0.5))), [0])[:] == [2/3]
@test freqresp(feedback(1.0, P2), ω)[:] ≈ 1 ./ (1 .+ P2_fr)

@test freqresp(feedback(0.5, delay(2.0)), ω) ≈ 0.5 ./ (1 .+ 0.5*exp.(-2im*ω))
@test freqresp(feedback(delay(2.0), 0.5), ω) ≈ exp.(-2im*ω) ./ (1 .+ 0.5*exp.(-2im*ω))

@test freqresp(feedback(P1, delay(1)), ω)[:] ≈ P1_fr ./ (1 .+ exp.(-im*ω) .* P1_fr) rtol=1e-15
@test freqresp(feedback(delay(1), P1), ω)[:] ≈ exp.(-im*ω) ./ (1 .+ exp.(-im*ω) .* P1_fr) rtol=1e-15 #FIXME: Answer is Inf, but should give error.. rtol=1e-15
@test freqresp(feedback(P1*delay(1), 1.0), ω)[:] ≈ P1_fr .* exp.(-im*ω) ./ (1 .+ exp.(-im*ω) .* P1_fr) rtol=1e-15
@test freqresp(feedback(1.0, P1*delay(1)), ω)[:] ≈ 1 ./ (1 .+ exp.(-im*ω) .* P1_fr) rtol=1e-15

@test freqresp(feedback(1.0, P2*0.5*(ss(1.0) + delay(2))), ω)[:] ≈ 1 ./(1 .+ P2_fr .* 0.5.*(1 .+ exp.(-2*im*ω)))


@test freqresp(1.0 + delay(2), ω)[:] ≈ 1 .+ exp.(-2im*ω) rtol=1e-15

G = ss(0.5) + 0.5*delay(2)# + 0.5*delay(3)
G_fr = 0.5 .+ 0.5*exp.(-2*im*ω)# + 0.5*exp.(-3*im*ω)
@test freqresp(G, ω)[:] ≈ G_fr rtol=1e-15

@test freqresp(feedback(1.0, P1*G), ω)[:] ≈ 1 ./(1 .+ P1_fr .* G_fr) rtol=1e-15
@test freqresp(feedback(P1, G), ω)[:] ≈ P1_fr ./(1 .+ P1_fr .* G_fr) rtol=1e-15

@test freqresp(feedback(P2, G), ω)[:] ≈ P2_fr ./(1 .+ P2_fr .* G_fr) rtol=1e-15

s = tf("s")

# Test alternative exp constructor for delays
d = exp(-2*s)
@test freqresp(d, [0, 1, 2]) ≈ [1, exp(-2im), exp(-4im)]

@test_throws ErrorException exp(-s^2 - 2*s)
@test_throws ErrorException exp(-2*s+1) # in principle ok, but not allowed anyway
@test_throws ErrorException exp([-2*s; -s])
@test_throws ErrorException exp(2*s) # Non-causal


Ω = [0.5, 1.0, 1.5, 2.0]
# Test for internal function delayd_ss
@test freqresp(ControlSystems.delayd_ss(1.0, 0.2), Ω)[:] ≈ exp.(-im*Ω) atol=1e-14
@test freqresp(ControlSystems.delayd_ss(3.2, 0.4), Ω)[:] ≈ exp.(-3.2*im*Ω) atol=1e-14
@test_throws ErrorException freqresp(ControlSystems.delayd_ss(3.2, 0.5), Ω)

# Simple tests for c2d of DelayLtiSystems
@test freqresp(c2d(feedback(ss(0,1,1,0), delay(1.5)), 0.5), Ω) ≈ [0.5/((z - 1) + 0.5*z^-3) for z in exp.(im*Ω*0.5)]
@test freqresp(c2d(feedback(delay(1.5), delay(1.0)), 0.5), Ω) ≈ [z^-3/(1 + z^-5) for z in exp.(im*Ω*0.5)]
@test freqresp(c2d(feedback(0.5, delay(1.5)), 0.5), Ω) ≈ [0.5/(1 + 0.5*z^-3) for z in exp.(im*Ω*0.5)]


# Random conversions
sys1 = DelayLtiSystem(1.0/s)
@test sys1.P.A == sys1.P.D == fill(0,1,1)
@test sys1.P.B*sys1.P.C == fill(1,1,1)
@test sys1.Tau == []

## Multiple Delays
G = ss(1.0) + delay(2) + delay(3)
G_fr = 1 .+ exp.(-2*im*ω) .+ exp.(-3*im*ω)
@test freqresp(G, ω)[:] ≈ G_fr rtol=1e-15
@test freqresp(feedback(1.0, G), ω)[:] ≈ 1 ./(1 .+ G_fr) # Somewhat pathological system though

sys = feedback(1.0, P1 * (delay(2) + delay(3)))
expected_sys_fr = 1.0 ./ (P1_fr .* (exp.(-2*im*ω) + exp.(-3*im*ω)) .+ 1)

@test freqresp(sys, ω)[:] ≈ expected_sys_fr rtol=1e-14

@test freqresp(feedback(1.0*P1, G*G*P2*P2), ω)[:] ≈ P1_fr ./(1 .+ (G_fr.^2).*P1_fr.*P2_fr.^2) rtol=1e-15# Somewhat pathological system though

s11 = feedback(ss(1/s), delay(1))
s12 = ss(1/s)
s21 = DelayLtiSystem(ss(1/(s+1)))
s22 = ss(1/(s+10))

s1 = [s11 s12]
s2 = [s21 s22]

f1 = [s1;s2]
f2 = [s11 s12;
      s21 s22]

# Test that different consatenations work
w = 10 .^ (-2:0.1:2)
@test freqresp(f1, w) ≈ freqresp(f2, w) rtol=1e-15

# This used to be a weird bug
@test freqresp(s11, w) ≈ freqresp(f2[1,1], w) rtol=1e-15


@test propertynames(delay(1.0)) == (:P, :Tau)


#FIXME: A lot more tests, including MIMO systems in particular

# Test step
println("Simulating first delay system:")
@time step(delay(1)*tf(1,[1.,1]))
@time step(delay(1)*tf(1,[1,1]))

@time y1, t1, x1 = step([s11;s12], 10)
@time @test y1[:,2] ≈ step(s12, t1)[1] rtol = 1e-14

t = 0.0:0.1:10
y2, t2, x2 = step(s1, t)
# TODO Figure out which is inexact here
@test y2[:,1,1:1] + y2[:,1,2:2] ≈ step(s11, t)[1] + step(s12, t)[1] rtol=1e-14

y3, t3, x3 = step([s11; s12], t)
@test y3[:,1,1] ≈ step(s11, t)[1] rtol = 1e-14
@test y3[:,2,1] ≈ step(s12, t)[1] rtol = 1e-14

y1, t1, x1 = step(DelayLtiSystem([1.0/s 2/s; 3/s 4/s]), t)
y2, t2, x2 = step([1.0/s 2/s; 3/s 4/s], t)
@test y1 ≈ y2 rtol=1e-14
@test size(x1,1) == length(t)
@test size(x1,3) == 2


##################### Test with known solution
K = 2.0
sys_known = feedback(delay(1)*K/s, 1)

ystep, t, _ = step(sys_known, 3)

function y_expected(t, K)
      if t < 1
            return 0
      elseif t < 2
            return K*(t-1)
      elseif t <= 3
            return K*(t-1)-1/2*K^2*(t-2)^2
      else
            throw(ArgumentError("Test not defined here"))
      end
end

@test ystep ≈ y_expected.(t, K) atol = 1e-12

function dy_expected(t, K)
      if t < 1
            return 0
      elseif t < 2
            return K
      elseif t <= 3
            return K - K^2*(t-2)
      else
            throw(ArgumentError("Test not defined here"))
      end
end

y_impulse, t, _ = impulse(sys_known, 3)

# TODO Better accuracy for impulse
@test y_impulse ≈ dy_expected.(t, K) rtol=1e-13
@test maximum(abs, y_impulse - dy_expected.(t, K)) < 1e-12

y_impulse, t, _ = impulse(sys_known, 3, alg=MethodOfSteps(Tsit5()))
# Two orders of magnitude better with BS3 in this case, which is default for impulse
@test y_impulse ≈ dy_expected.(t, K) rtol=1e-5
@test maximum(abs, y_impulse - dy_expected.(t, K)) < 1e-5

## Test delay with D22 term
t = 0:0.01:4

sys = delay(1)

y, t, x = step(sys, t)
@test y[:] ≈ [zeros(100); ones(301)] atol = 1e-12
@test size(x) == (401,0)

sys = delay(1)*delay(1)*1/s

y, t, x = step(sys, t)

y_sol = [zeros(200);0:0.01:2]

@test maximum(abs,y-y_sol) < 1e-13
@test maximum(abs,x-collect(0:0.01:4)) < 1e-15

# TODO For some reason really bad accuracy here
# Looks like a lag in time
sys = 1/s*delay(1)*delay(1)

y, t, x = step(sys, t)
@test maximum(abs,y-y_sol) < 1e-5
@test maximum(abs,x-y_sol) < 1e-5

t = 0:0.001:0.1
y, t, x = step(sys, t)
@test length(y) == length(t)

##  Test of basic pade functionality

Ω = [0, 0.5, 1, 2, 5]
@test freqresp(pade(1, 1), Ω) == freqresp(tf([-1/2, 1], [1/2, 1]), Ω)
@test freqresp(pade(1, 2), Ω) ≈ freqresp(tf([1/12, -1/2, 1], [1/12, 1/2, 1]), Ω)

for (n, tol)=enumerate([0.05; 1e-3; 1e-5; 1e-7; 1e-11; 1e-14*ones(5)])
    G = pade(0.8, n)

    @test isstable(G)
    @test evalfr(G, 0)[1] ≈ 1
    @test abs(evalfr(G, 2im)[1]) ≈ 1
    @test evalfr(G, 1im)[1] ≈ exp(-0.8im) atol=tol
end


## Test pade applied to DelayLtiSystem

@test freqresp(pade(delay(0.5), 2), Ω) ≈ freqresp(pade(0.5, 2), Ω)

P = delay(1.0) * DemoSystems.lag(T=1)
@test freqresp(pade(feedback(1,P), 2), Ω) == freqresp(feedback(1, pade(P,2)), Ω)


Ω = [0, 0.1, 0.2]
P_wb = DemoSystems.woodberry()

@test freqresp(pade(P_wb, 2), Ω) ≈ freqresp(P_wb, Ω) atol=0.02
@test freqresp(pade(P_wb, 3), Ω) ≈ freqresp(P_wb, Ω) atol=5e-4

@test freqresp(pade(feedback(eye_(2), P_wb), 3), Ω) ≈ freqresp(feedback(eye_(2), P_wb), Ω) atol=1e-4

end
