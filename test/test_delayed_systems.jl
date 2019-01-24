ω = 0.0:8

# broken: typeof(promote(delay(0.2), ss(1))[1]) == DelayLtiSystem{Float64}

@test typeof(promote(delay(0.2), ss(1.0 + im))[1]) == DelayLtiSystem{Complex{Float64}}

@test freqresp(delay(1), ω)[:] ≈ exp.(-im*ω) rtol=1e-15
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

# Addition
@test_broken freqresp(1 + delay(1), ω)[:] ≈ 1 .+ exp.(-im*ω) # FIXME; Add suitable conversion from Int64 to DelayLtiSystem
@test freqresp(P1 + delay(1), ω)[:] ≈ P1_fr .+ exp.(-im*ω)

#FIXME: the following gives a crash.. freqresp(P1 - delay(1), ω)[:] ≈ P1_fr .+ exp.(-im*ω)


# Multiplication
@test freqresp(P1 * delay(1), ω)[:] ≈ P1_fr .* exp.(-im*ω) rtol=1e-15
@test freqresp(delay(1) * P1, ω)[:] ≈ P1_fr .* exp.(-im*ω) rtol=1e-15

@test freqresp(P2 * delay(1), ω)[:] ≈ P2_fr .* exp.(-im*ω) rtol=1e-15
@test freqresp(delay(1) * P2, ω)[:] ≈ P2_fr .* exp.(-im*ω) rtol=1e-15


# Feedback
# The first tests don't include delays, but the linear system is of DelayLtiForm type
# (very simple system so easy to troubleshoot)
@test freqresp(feedback(1.0, P1), ω)[:] ≈ 1 ./ (1 .+  P1_fr) rtol=1e-15
@test freqresp(feedback(P1, 1.0), ω)[:] ≈ P1_fr ./ (1 .+  P1_fr) rtol=1e-15
@test freqresp(feedback(P1, P1), ω)[:] ≈ P1_fr ./ (1 .+ P1_fr .* P1_fr) rtol=1e-15

@test freqresp(feedback(1.0, DelayLtiSystem(ss(0.5))), [0])[:] == [2/3]
@test freqresp(feedback(1.0, P2), ω)[:] ≈ 1 ./ (1 .+ P2_fr)


freqresp(feedback(1.0, ss(0.5) + 0.5*delay(2) + 0.5*delay(3)), ω)[:]
freqresp(feedback(1.0, G), ω)[:]


@test freqresp(feedback(0.5, delay(2.0)), ω) ≈ 0.5 ./ (1 .+ 0.5*exp.(-2im*ω))
@test freqresp(feedback(delay(2.0), 0.5), ω) ≈ exp.(-2im*ω) ./ (1 .+ 0.5*exp.(-2im*ω))

@test freqresp(feedback(P1, delay(1)), ω)[:] ≈ P1_fr ./ (1 .+ exp.(-im*ω) .* P1_fr) rtol=1e-15
@test freqresp(feedback(delay(1), P1), ω)[:] ≈ exp.(-im*ω) ./ (1 .+ exp.(-im*ω) .* P1_fr) #FIXME: Answer is Inf, but should give error.. rtol=1e-15
@test freqresp(feedback(P1*delay(1), 1.0), ω)[:] ≈ P1_fr .* exp.(-im*ω) ./ (1 .+ exp.(-im*ω) .* P1_fr) rtol=1e-15
@test freqresp(feedback(1.0, P1*delay(1)), ω)[:] ≈ 1 ./ (1 .+ exp.(-im*ω) .* P1_fr) rtol=1e-15

@test freqresp(feedback(1.0, P2*0.5*(ss(1.0) + delay(2))), ω)[:] ≈ 1 ./(1 .+ P2_fr .* 0.5.*(1 .+ exp.(-2*im*ω)))


@test_broken freqresp(1.0 + delay(2), ω)[:] # TODO: Add addition for DelayLtiSystem. how to do this in a convenient ammner?

G = ss(0.5) + 0.5*delay(2)# + 0.5*delay(3)
G_fr = 0.5 .+ 0.5*exp.(-2*im*ω)# + 0.5*exp.(-3*im*ω)
@test freqresp(G, ω)[:] ≈ G_fr rtol=1e-15

@test freqresp(feedback(1.0, P1*G), ω)[:] ≈ 1 ./(1 .+ P1_fr .* G_fr) rtol=1e-15
@test freqresp(feedback(P1, G), ω)[:] ≈ P1_fr ./(1 .+ P1_fr .* G_fr) rtol=1e-15

@test freqresp(feedback(P2, G), ω)[:] ≈ P2_fr ./(1 .+ P2_fr .* G_fr) rtol=1e-15

sys = feedback(1.0, P1 * (delay(2) + delay(3)))
expected_sys_fr = 1.0 ./ (P1_fr .* (exp.(-2*im*ω) + exp.(-3*im*ω)) .+ 1)
@test freqresp(sys, ω)[:] ≈ expected_sys_fr rtol=1e-14

## Multiple Delays
G = ss(1.0) + delay(2) + delay(3)
G_fr = 1 .+ exp.(-2*im*ω) .+ exp.(-3*im*ω)
@test freqresp(G, ω)[:] ≈ G_fr rtol=1e-15
@test freqresp(feedback(1.0, G), ω)[:] ≈ 1 ./(1 .+ G_fr) # Somewhat pathological system though


#FIXME: A lot more tests, including MIMO systems in particular
