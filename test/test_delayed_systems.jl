@testset "test_delay_system" begin
ω = 0.0:8

# broken: typeof(promote(delay(0.2), ss(1))[1]) == DelayLtiSystem{Float64}

@test typeof(promote(delay(0.2), ss(1.0 + im))[1]) == DelayLtiSystem{Complex{Float64}}

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
@test freqresp(1 + delay(1), ω)[:] ≈ 1 .+ exp.(-im*ω) # FIXME; Add suitable conversion from Int64 to DelayLtiSystem
@test freqresp(P1 + delay(1), ω)[:] ≈ P1_fr .+ exp.(-im*ω)

# Substraction
@test freqresp(P1 - delay(1), ω)[:] ≈ P1_fr .- exp.(-im*ω)


## Multiplication
@test freqresp(P1 * delay(1), ω)[:] ≈ P1_fr .* exp.(-im*ω) rtol=1e-15
@test freqresp(delay(1) * P1, ω)[:] ≈ P1_fr .* exp.(-im*ω) rtol=1e-15

@test freqresp(P2 * delay(1), ω)[:] ≈ P2_fr .* exp.(-im*ω) rtol=1e-15
@test freqresp(delay(1) * P2, ω)[:] ≈ P2_fr .* exp.(-im*ω) rtol=1e-15


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


#FIXME: A lot more tests, including MIMO systems in particular

# Test step
y1, t1, x1 = step([s11;s12], 10)
@test y1[:,2] ≈ step(s12, t1)[1] rtol = 1e-14

t = 0.0:0.1:10
y2, t2, x2 = step(s1, t)
# TODO Figure out which is inexact here
@test y2[:,1,1:1] + y2[:,1,2:2] ≈ step(s11, t)[1] + step(s12, t)[1] rtol=1e-5

y3, t3, x3 = step([s11; s12], t)
@test y3[:,1,1] ≈ step(s11, t)[1] rtol = 1e-4
@test y3[:,2,1] ≈ step(s12, t)[1] rtol = 1e-14

y1, t1, x1 = step(DelayLtiSystem([1.0/s 2/s; 3/s 4/s]), t)
y2, t2, x2 = step([1.0/s 2/s; 3/s 4/s], t)
@test y1 ≈ y2 rtol=1e-15
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

@test ystep ≈ y_expected.(t, K) atol = 1e-13

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
@test y_impulse ≈ dy_expected.(t, K) rtol=1e-2
@test maximum(abs, y_impulse - dy_expected.(t, K)) < 1e-2

@time [s11; s12]

end
