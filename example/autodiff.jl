using ControlSystems, OrdinaryDiffEq
p0          = [1., 1., 1.]             # optimeringsvariabel
K(kp,ki,kd) = pid(kp=kp, ki=ki, kd=kd) # K regulator, s Laplacevariabel (alltihop är i kontinuerlig tid)
K(p)        = K(p...)

ζ  = 0.3
ω  = 1.; ω² = ω^2
P  = tf(ω²,[1, 2ζ*ω, ω²])^2 # LTI-process, tillhandahålls av användaren
# P = c2d(P,h)[1]
Ω  = logspace(-3,3,100) # vektor med vinkelfrekvenser för utvärdering av bivillkor, tillhandahålls av användaren
h  = 0.5
Tf = 30.
t  = 0:h:Tf
 # Ms och Mt är positiva skalärer, tillhandahålls av användaren
Ms = 1.2
Mt = 2

p  = p0
# Diverse system
function costfun(p)
    C = K(p)
    # C = tf(1,h)
    L = P*C
    S = feedback(tf(1.),L)
    T = 1-S
    PS = ss(P*S)
    C = ss(C)
    s = Simulator(PS, (t,x) -> [1])
    sol = solve(s, zeros(PS.nx), (0.,Tf), Tsit5())
    y = PS.C*sol(t)
    Sw = vec(bode(S,Ω)[1])
    Tw = vec(bode(T,Ω)[1])
    J = sum(abs2,1-y) # kostnaden som ska minimeras är L1-norm av PS stegsvar
    # y,_,_ = step(PS, t)
    J += 100*sum((Sw.>Ms))
    J += 100*sum((Tw.>Mt))
    J
end
costfun(p0)


res = bboptimize(costfun, NumDimensions = 3, MaxTime = 30.0)
res = compare_optimizers(costfun; SearchRange = (0.,10.), NumDimensions = 3, MaxTime = 20.0)
p = best_candidate(res)

using ForwardDiff
ForwardDiff.gradient(costfun, p)

using ReverseDiff
ReverseDiff.gradient(costfun, p)

#
# p >= 0                                          # elementvis
# cS(p) = abs(S(p) - Ms <= 0              # vektorvärt, med ett element för varje element i w
# cT(p) = abs(T(p) - Mt <= 0              # samma
