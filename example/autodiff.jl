using ControlSystems, OrdinaryDiffEq
p0          = [0.1, 0.1, 0.1]             # optimeringsvariabel
K(kp,ki,kd) = pid(kp=kp, ki=ki, kd=kd)
K(p)        = K(p...)

ζ  = 0.3
ω  = 1.; ω² = ω^2
P  = tf(ω²,[1, 2ζ*ω, ω²])^2 # LTI-process, tillhandahålls av användaren

Ω  = logspace(-1,1,100) # vektor med vinkelfrekvenser för utvärdering av bivillkor, tillhandahålls av användaren
h  = 0.1
Tf = 30.
t  = 0:h:Tf-h
 # Ms och Mt är positiva skalärer, tillhandahålls av användaren
Ms = 1.5
Mt = 1.8

p  = copy(p0)

function costfun(p)
    C     = K(p)
    L     = P*C
    S     = feedback(tf(1.),L)
    T     = 1-S
    PS    = ss(P*S)
    s     = Simulator(PS, (t,x) -> [1])
    x0    = zeros(PS.nx)
    ty    = eltype(p)
    tspan = (ty(0.),ty(Tf+5))
    sol   = solve(s, ty.(x0), tspan, Tsit5())
    y     = PS.C*sol(t)
    Sw    = vec(bode(S,Ω)[1])
    Tw    = vec(bode(T,Ω)[1])
    J     = sum(abs,1-y) # kostnaden som ska minimeras är L1-norm av PS stegsvar
    J -= sum(x->log(0.1*max(x,1e-5)),Ms-Sw)  #(sum((Sw.>Ms))+1)
    J -= sum(x->log(0.1*max(x,1e-5)),Mt-Tw) #(sum((Tw.>Mt))+1)
    J
end
costfun(p0)


using Optim
od = OnceDifferentiable(costfun, p; autodiff = :forward)
td = TwiceDifferentiable(costfun, p; autodiff = :forward)
@time res3 = optimize(od, p, BFGS(), Optim.Options(show_trace=true))
@time res4 = optimize(td, p, Newton(), Optim.Options(show_trace=true))
p = Optim.minimizer(res3)

using BlackBoxOptim
res1 = compare_optimizers(costfun; SearchRange = (0.,1.), NumDimensions = 3, MaxTime = 20.0)
res2 = bboptimize(costfun, NumDimensions = 3, MaxTime = 30.0)
p = best_candidate(res2)
# optimize(costfun, p)


#
# p >= 0                                          # elementvis
# cS(p) = abs(S(p) - Ms <= 0              # vektorvärt, med ett element för varje element i w
# cT(p) = abs(T(p) - Mt <= 0              # samma


using ReverseDiff
ReverseDiff.gradient(costfun, p)
