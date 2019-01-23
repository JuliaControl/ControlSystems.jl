ω = 0.0:8

typeof(promote(delay(0.2), ss(1))[1]) == DelayLtiSystem{Float64}
typeof(promote(delay(0.2), ss(1.0 + im))[1]) == DelayLtiSystem{Float64}

freqresp(delay(1), ω)[:] ≈ exp.(-im*ω)
freqresp(delay(2.5), ω)[:] ≈ exp.(-2.5*im*ω)

P1 = DelayLtiSystem(ss(-1.0, 1, 1, 0))
P1_fr = 1 ./ (im*ω .+ 1)

P2 = DelayLtiSystem(ss(-2.0, -1, 1, 1)) # (s+1)/(s+2)
P2_fr = (im*ω .+ 1) ./ (im*ω .+ 2)


# Addition
freqresp(P1 + delay(1), ω)[:] ≈ P1_fr .+ exp.(-im*ω)
#FIXME: the following gives a crash.. freqresp(P1 - delay(1), ω)[:] ≈ P1_fr .+ exp.(-im*ω)


# Multiplication
freqresp(P1 * delay(1), ω)[:] ≈ P1_fr .* exp.(-im*ω)
freqresp(delay(1) * P1, ω)[:] ≈ P1_fr .* exp.(-im*ω)

freqresp(P2 * delay(1), ω)[:] ≈ P2_fr .* exp.(-im*ω)
freqresp(delay(1) * P2, ω)[:] ≈ P2_fr .* exp.(-im*ω)


# Feedback
# The first tests don't include delays, but the linear system is of DelayLtiForm type
# (very simple system so easy to troubleshoot)
freqresp(feedback(1.0, P1), ω)[:] ≈ 1 ./ (1 .+  P1_fr)
freqresp(feedback(P1, 1.0), ω)[:] ≈ P1_fr ./ (1 .+  P1_fr)
freqresp(feedback(P1, P1), ω)[:] ≈ P1_fr ./ (1 .+ P1_fr .* P1_fr)

freqresp(feedback(P1, delay(1)), ω)[:] ≈ P1_fr ./ (1 .+ exp.(-im*ω) .* P1_fr)
freqresp(feedback(delay(1), P1), ω)[:] ≈ exp.(-im*ω) ./ (1 .+ exp.(-im*ω) .* P1_fr) #FIXME: Answer is Inf, but should give error..
freqresp(feedback(P1*delay(1), 1.0), ω)[:] ≈ P1_fr .* exp.(-im*ω) ./ (1 .+ exp.(-im*ω) .* P1_fr)
freqresp(feedback(1.0, P1*delay(1)), ω)[:] ≈ 1 ./ (1 .+ exp.(-im*ω) .* P1_fr)

#FIXME: A lot more tests, including MIMO systems in particular
