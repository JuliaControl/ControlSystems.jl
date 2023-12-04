using ControlSystemsBase

@testset "test_delay_system" begin

# For simplicity, equality of DelayLtiSystems are tested over a finite set of frequencies
ω = 0.0:8


@test typeof(promote(delay(0.2), ss(1))[1]) == DelayLtiSystem{Float64,Float64}

@test typeof(promote(delay(0.2), ss(1.0 + im))[1]) == DelayLtiSystem{Complex{Float64}, Float64}

@test delay(1, 1) isa ControlSystemsBase.StateSpace{<:Discrete}
@test_throws ErrorException delay(1, 0.3)
@test delay(1, 1) == ss(zeros(1,1), ones(1,1), ones(1,1), zeros(1,1), 1)

@test sprint(show, ss(1,1,1,1)*delay(1.0)) == "DelayLtiSystem{Float64, Float64}\n\nP: StateSpace{Continuous, Float64}\nA = \n 1.0\nB = \n 0.0  1.0\nC = \n 1.0\n 0.0\nD = \n 0.0  1.0\n 1.0  0.0\n\nContinuous-time state-space model\n\nDelays: [1.0]"

# Extremely basic tests
@test freqresp(delay(1), ω) ≈ reshape(exp.(-im*ω), 1, 1, length(ω)) rtol=1e-15
@test freqresp(delay(2.5), ω)[:] ≈ exp.(-2.5im*ω) rtol=1e-15
@test freqresp(3.5*delay(2.5), ω)[:] ≈ 3.5*exp.(-2.5im*ω) rtol=1e-15
@test freqresp(delay(2.5)*1.5, ω)[:] ≈ exp.(-2.5im*ω)*1.5 rtol=1e-15

# Addition of constant
@test evalfr(1 + delay(1.0), 0)[] ≈ 2
@test evalfr(1 - delay(1.0), 0)[] ≈ 0
@test evalfr([2 -delay(1.0)], 0) ≈ [2 -1]

# Strictly proper system
P1 = DelayLtiSystem(ss(-1.0, 1, 1, 0))
P1_fr = 1 ./ (im*ω .+ 1)
@test freqresp(P1, ω)[:] ≈ P1_fr

# Not stritcly proper system
P2 = DelayLtiSystem(ss(-2.0, -1, 1, 1)) # (s+1)/(s+2)
P2_fr = (im*ω .+ 1) ./ (im*ω .+ 2)
@test freqresp(P2, ω)[:] ≈ P2_fr rtol=1e-15


## Addition
@test freqresp(1 + delay(1), ω)[:] ≈ 1 .+ exp.(-im*ω)
@test freqresp(P1 + delay(1), ω)[:] ≈ P1_fr .+ exp.(-im*ω)

# Subtraction
@test freqresp(P1 - delay(1), ω)[:] ≈ P1_fr .- exp.(-im*ω)

## Multiplication by scalar
@test freqresp(2*delay(1), ω)[:] ≈ 2*exp.(-im*ω)
@test freqresp(2*(P1*delay(1)), ω)[:] ≈ 2*P1_fr.*exp.(-im*ω)

## Multiplication
@test freqresp(P1 * delay(1), ω)[:] ≈ P1_fr .* exp.(-im*ω) rtol=1e-15
@test freqresp(delay(1) * P1, ω)[:] ≈ P1_fr .* exp.(-im*ω) rtol=1e-15

@test freqresp(P2 * delay(1), ω)[:] ≈ P2_fr .* exp.(-im*ω) rtol=1e-15
@test freqresp(delay(1) * P2, ω)[:] ≈ P2_fr .* exp.(-im*ω) rtol=1e-15

## Division / feedback
@test freqresp(1/(1+P1), ω) ≈ freqresp(feedback(I(size(P1, 1)), P1), ω) rtol=1e-15

## Zero
@test zero(typeof(P1)) == ss(0.0)
@test zero(P1) == ss(0.0)

## append
P12 = append(P1, P2)
G12 = [P1 tf(0); tf(0) P2]

F = freqresp(P12, ω)

@test freqrespv(P1, ω) ≈ F[1,1,:]
@test freqrespv(P2, ω) ≈ F[2,2,:]
@test all(iszero, F[1,2,:])
@test all(iszero, F[2,1,:])

F = freqresp(G12, ω)
@test freqrespv(P1, ω) ≈ F[1,1,:]
@test freqrespv(P2, ω) ≈ F[2,2,:]
@test all(iszero, F[1,2,:])
@test all(iszero, F[2,1,:])

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

@test freqrespv(feedback(0.5, delay(2.0)), ω) ≈ 0.5 ./ (1 .+ 0.5*exp.(-2im*ω))
@test freqrespv(feedback(delay(2.0), 0.5), ω) ≈ exp.(-2im*ω) ./ (1 .+ 0.5*exp.(-2im*ω))

@test freqrespv(feedback(P1, delay(1)), ω) ≈ P1_fr ./ (1 .+ exp.(-im*ω) .* P1_fr) rtol=1e-15
@test freqrespv(feedback(delay(1), P1), ω) ≈ exp.(-im*ω) ./ (1 .+ exp.(-im*ω) .* P1_fr) rtol=1e-15 #FIXME: Answer is Inf, but should give error.. rtol=1e-15
@test freqrespv(feedback(P1*delay(1), 1.0), ω) ≈ P1_fr .* exp.(-im*ω) ./ (1 .+ exp.(-im*ω) .* P1_fr) rtol=1e-15
@test freqrespv(feedback(1.0, P1*delay(1)), ω) ≈ 1 ./ (1 .+ exp.(-im*ω) .* P1_fr) rtol=1e-15

@test freqrespv(feedback(1.0, P2*0.5*(ss(1.0) + delay(2))), ω) ≈ 1 ./(1 .+ P2_fr .* 0.5.*(1 .+ exp.(-2*im*ω)))
@test freqrespv(feedback(1.0, 0.5*P2*(ss(1.0) + delay(2))), ω) ≈ 1 ./(1 .+ P2_fr .* 0.5.*(1 .+ exp.(-2*im*ω)))


@test freqrespv(1.0 + delay(2), ω) ≈ 1 .+ exp.(-2im*ω) rtol=1e-15

G = ss(0.5) + 0.5*delay(2)# + 0.5*delay(3)
G_fr = 0.5 .+ 0.5*exp.(-2*im*ω)# + 0.5*exp.(-3*im*ω)
@test freqrespv(G, ω) ≈ G_fr rtol=1e-15

@test freqrespv(feedback(1.0, P1*G), ω) ≈ 1 ./(1 .+ P1_fr .* G_fr) rtol=1e-15
@test freqrespv(feedback(P1, G), ω) ≈ P1_fr ./(1 .+ P1_fr .* G_fr) rtol=1e-15

@test freqrespv(feedback(P2, G), ω) ≈ P2_fr ./(1 .+ P2_fr .* G_fr) rtol=1e-15

s = tf("s")

# Test alternative exp constructor for delays
d = exp(-2*s)
@test freqrespv(d, [0, 1, 2]) ≈ [1, exp(-2im), exp(-4im)]

@test_throws ErrorException exp(-s^2 - 2*s)
@test_throws ErrorException exp(-2*s+1) # in principle ok, but not allowed anyway
@test_throws ErrorException exp([-2*s; -s])
@test_throws ErrorException exp(2*s) # Non-causal


Ω = [0.5, 1.0, 1.5, 2.0]
# Test for internal function delayd_ss
@test freqresp(ControlSystemsBase.delayd_ss(1.0, 0.2), Ω)[:] ≈ exp.(-im*Ω) atol=1e-14
@test freqresp(ControlSystemsBase.delayd_ss(3.2, 0.4), Ω)[:] ≈ exp.(-3.2*im*Ω) atol=1e-14
@test_throws ErrorException ControlSystemsBase.delayd_ss(3.2, 0.5)

# Simple tests for c2d of DelayLtiSystems
@test freqrespv(c2d(feedback(ss(0,1,1,0), delay(1.5)), 0.5), Ω) ≈ [0.5/((z - 1) + 0.5*z^-3) for z in exp.(im*Ω*0.5)]
@test freqrespv(c2d(feedback(delay(1.5), delay(1.0)), 0.5), Ω) ≈ [z^-3/(1 + z^-5) for z in exp.(im*Ω*0.5)]
@test freqrespv(c2d(feedback(0.5, delay(1.5)), 0.5), Ω) ≈ [0.5/(1 + 0.5*z^-3) for z in exp.(im*Ω*0.5)]


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

# Test that different concatenations work
w = 10 .^ (-2:0.1:2)
@test freqresp(f1, w) ≈ freqresp(f2, w) rtol=1e-15

# This used to be a weird bug
@test freqresp(s11, w) ≈ freqresp(f2[1,1], w) rtol=1e-15


@test propertynames(delay(1.0)) == (:P, :Tau, :nu, :ny)

@test_throws SingularException 1/f2
@test_throws SingularException randn(2,2)/f2
@test_throws SingularException f2/f2
@test 1/(I(2)+f2) == feedback(I(2), f2)

#FIXME: A lot more tests, including MIMO systems in particular

# Test step
@test_throws MethodError step(delay(1)*tf(1,[1.,1]))


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


# test automatic frequency selection
mag, phase, w = bode(DemoSystems.lag()*delay(1))
@test w[1] <= 0.05
@test w[end] >= 5


end
