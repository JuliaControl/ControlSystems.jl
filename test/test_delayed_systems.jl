ω = 0.0:8

freqresp(delay(1), ω)[:] ≈ exp.(-im*ω)
freqresp(delay(2.5), ω)[:] ≈ exp.(-2.5*im*ω)

P1 = DelayLtiSystem(ss(-1.0, 1, 1, 0))
P1_fr ≈ 1 ./ (im*ω .+ 1)

P2 = DelayLtiSystem(ss(-2.0, -1, 1, 1)) # (s+1)/(s+2)
P2_fr = (im*ω .+ 1) ./ (im*ω .+ 2)


freqresp(P1 * delay(1), ω)[:] ≈ P1_fr .* exp.(-im*ω)
freqresp(delay(1) * P1, ω)[:] ≈ P1_fr .* exp.(-im*ω)

freqresp(P2 * delay(1), ω)[:] ≈ P2_fr .* exp.(-im*ω)
freqresp(delay(1) * P2, ω)[:] ≈ P2_fr .* exp.(-im*ω)


# feedback
freqresp(feedback(P1, delay(1)), ω)[:] ≈ P1_fr ./ (1 .+ exp.(-im*ω) .* P1_fr)
freqresp(feedback(delay(1), P1), ω)[:] ≈ exp.(-im*ω) ./ (1 .+ exp.(-im*ω) .* P1_fr)

G_fr = freqresp(feedback(P1*delay(1), DelayLtiSystem(ss(1.0))), ω)[:]
G_fr ≈ P1_fr .* exp.(-im*ω) ./ (1 .+ exp.(-im*ω) .* P1_fr)


# Test DelayLtiSystem objects without delays
freqresp(feedback(P1, P1), ω)[:] ≈ P1_fr ./ (1 .+ P1_fr .* P1_fr)


# FIXME: Automatic conversion to DelayLtiSystem


nvert(DelayLtiSystem, P1)

promote_rule(typeof(P1), typeof(P2))
x1, x2 = promote(P1, P2)


# What ordering for freqresp ?!
freqresp(P.P, ω)

G_fr = freqresp(G, ω)[:]
