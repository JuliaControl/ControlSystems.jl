module TestAbstractTF
using CustomTest
using ControlSystems

# CONTINUOUS
@test C_011 = tf("(s+2)")
@test C_111 = tf("(s+1)/(s+2)")
@test C_111*C_011 == tf("(s+2)*((s+1)/(s+2))")

@test C_011(1im) == 2+1im
@test C_111([1im, 2im]) == [0.6+0.2im, 0.75 + 0.25im]
@test (C_111*C_011)(im) == 1.0+1im

bode(C_011)
end
