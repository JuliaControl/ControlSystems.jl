module TestGeneralizedTF
using CustomTest
using ControlSystems

# CONTINUOUS
C_011 = tfg("(s+2)")
C_111 = tfg("(s+1)/(s+2)")
@test C_011*C_111 == tfg("(s+2)*((s+1)/(s+2))")

#We might want to make evalfr scalar
@test C_011(1im) == reshape([2+1im;],1,1)
@test (C_111*C_011)(im) == reshape([1.0+1im],1,1)

@test tf(C_111) == tf([1,1],[1,2])
@test zpk(C_011*C_111) == zpk([-1,-2],[-2],1)

@test bode(C_111*C_011, logspace(-1,1)) == bode(tfg("(s+2)*((s+1)/(s+2))"), logspace(-1,1))

# Test numpoly, numvec, denpoly, denvec for SisoZpk

# Test deprecation (is not an error)
#@test_err num(C_111.matrix[1,1])
#@test_err den(C_111.matrix[1,1])

@test_err numvec(C_111.matrix[1,1])
@test_err denvec(C_111.matrix[1,1])

@test_err numvec(C_111)
@test_err denvec(C_111)

@test_err numpoly(C_111)
@test_err denpoly(C_111)

end
