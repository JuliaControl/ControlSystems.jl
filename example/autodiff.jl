using ControlSystems, OrdinaryDiffEq
p0          = [0.01, 0.01, 0.01]             # optimeringsvariabel
K(kp,ki,kd) = pid(kp=kp, ki=ki, kd=kd)
K(p)        = K(p...)

ζ  = 0.1
ω  = 1.; ω² = ω^2
P  = tf(ω²,[1, 2ζ*ω, ω²])^2 # LTI-process, tillhandahålls av användaren

Ω  = logspace(-1,2,150) # vektor med vinkelfrekvenser för utvärdering av bivillkor, tillhandahålls av användaren
h  = 0.1
Tf = 60.
t  = 0:h:Tf-h
# Ms och Mt är positiva skalärer, tillhandahålls av användaren
Ms = 1.7
Mt = 1.7

p  = copy(p0)

function costfun(p)
    C     = K(p)
    S     = 1/(1+P*C)
    PS    = ss(P*S)
    T     = tf(1.) - S
    s     = Simulator(PS, (t,x) -> [1])
    x0    = zeros(PS.nx)
    ty    = eltype(p)
    tspan = (ty(0.),ty(Tf))
    sol   = solve(s, ty.(x0), tspan, Tsit5())
    y     = PS.C*sol(t)
    Sw    = vec(bode(S,Ω)[1])
    Tw    = vec(bode(T,Ω)[1])
    J0     = mean(abs,y) # kostnaden som ska minimeras är L1-norm av PS stegsvar
    J = J0
    # J *=  (sum((Sw.>Ms))+1)
    # J *=  (sum((Tw.>Mt))+1)
    J -= 0.1*sum(x->log(0.11*max(x,1e-5)),Ms-Sw)  #(sum((Sw.>Ms))+1)
    J -= 0.1*sum(x->log(0.11*max(x,1e-5)),Mt-Tw) #(sum((Tw.>Mt))+1)
    J
end
costfun(p0)


using Optim
od = OnceDifferentiable(costfun, p; autodiff = :forward)
@time res3 = optimize(od, p, BFGS(), Optim.Options(show_trace=true))
evalsol(res3)

td = TwiceDifferentiable(costfun, p; autodiff = :forward)
@time res4 = optimize(td, p, Newton(), Optim.Options(show_trace=true))

using BlackBoxOptim
res1 = compare_optimizers(costfun; SearchRange = (0.,2.), NumDimensions = length(p0), MaxTime = 20.0)
res2 = bboptimize(costfun, NumDimensions = length(p0), MaxTime = 15.0)
evalsol(res2)
p = best_candidate(res2)
# optimize(costfun, p)


#
# p >= 0                                          # elementvis
# cS(p) = abs(S(p) - Ms <= 0              # vektorvärt, med ett element för varje element i w
# cT(p) = abs(T(p) - Mt <= 0              # samma


using ReverseDiff
ReverseDiff.gradient(costfun, p)


evalsol(res::Optim.MultivariateOptimizationResults) = evalsol(Optim.minimizer(res))
evalsol(res::BlackBoxOptim.OptimizationResults) = evalsol(best_candidate(res))
function evalsol(p::Vector)
    C     = K(p)
    S     = 1/(1+P*C)
    PS    = ss(P*S)
    T     = tf(1.) - S
    s     = Simulator(PS, (t,x) -> [1])
    x0    = zeros(PS.nx)
    ty    = eltype(p)
    tspan = (ty(0.),ty(Tf))
    sol   = solve(s, ty.(x0), tspan, Tsit5())
    y     = PS.C*sol(t)
    Sw    = vec(bode(S,Ω)[1])
    Tw    = vec(bode(T,Ω)[1])

    plot(t,y', layout=2, show=false)
    plot!(Ω, [Sw Tw] , lab=["Sw" "Tw"], subplot=2, xscale=:log10, yscale=:log10, show=false)
    plot!([Ω[1],Ω[end]], [Ms,Ms], c = :black, l=:dash, subplot=2, show=false)
    plot!([Ω[1],Ω[end]], [Mt,Mt], c = :purple, l=:dash, subplot=2)
    gui()
end
