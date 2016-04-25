module TestGeneralizedTF
using CustomTest
using ControlSystems

# CONTINUOUS
C_011 = tfg("(s+2)")
C_111 = tfg("(s+1)/(s+2)")
@test C_011*C_111 == tfg("(s+2)*((s+1)/(s+2))")

#We might want to make evalfr scalar
@test C_011(1im) == 2+1im
@test (C_111*C_011)(im) == 1.0+1im

@test tf(C_111) == tf([1,1],[1,2])
@test zpk(C_011*C_111) == zpk([-1,-2],[-2],1)

@test bode(C_111*C_011, logspace(-1,1)) == bode(tfg("(s+2)*((s+1)/(s+2))"), logspace(-1,1))
end
