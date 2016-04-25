module TestConversion
using CustomTest
using ControlSystems

f11 = tf(ss([-1 0;1 1],[1;0],[1 1],0))
f12 = tf([1,0],[1,0,-1])

#Direct term was previously missing
f21 = tf(ss([-1 0;1 1],[1;0],[1 1],1))
f22 = tf([1,1,-1],[1,0,-1])

#Test ss2tf
@test_approx_eq_eps f11 f12 1e-15
@test_approx_eq_eps f21 f22 1e-15

#Test tf2ss
#ss is not unique so go back and forth
@test_approx_eq_eps f11 tf(ss(f12)) 1e-15
@test_approx_eq_eps f21 tf(ss(f22)) 1e-15

end
