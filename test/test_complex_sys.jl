module TestComplexSys
using CustomTest
using Base.Test
using ControlSystems

# Just some test for early prototyping of statespace systems of
# different numerical types

# Simple first order system
G1 = ss(-1-1im, 1, 1, 0)
w = logspace(-1,1)
G1_fr,_ = freqresp(G1, w)

# Test freqeuncy reposne
@test_approx_eq G1_fr 1 ./ (1im*w + 1 + 1im)

# Step response
t = 0:0.1:3
y, _, _ = step(G1, t)
@test_approx_eq y (1 - exp((-1-im)*t))/(1 + im)

# Poles of the system
@test pole(G1) == [-1-1im]


@test ss(1,1,1,0.5im) + 0.5im == ss(1,1,1,1.0im)


# Just try to add a transfer function with a complex
# state-space model to see that nothing crashes
s = tf("s")
ss(1.0 + 1im,1,1,1.0) + 1/(s+1)


# Add test for printing complex system


end
