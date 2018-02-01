using ControlSystems, OrdinaryDiffEq
p0          = [0.01, 0.01, 0.01] # Initial guess
K(kp,ki,kd) = pid(kp=kp, ki=ki, kd=kd)
K(p)        = K(p...)

# Define process model
ζ  = 0.1
ω  = 1.; ω² = ω^2
P  = tf(ω²,[1, 2ζ*ω, ω²])^2

Ω  = logspace(-1,2,150)  # Frequency vector to eval constraints
h  = 0.1 # Sample time for time-domain evaluation
Tf = 60.  # Time horizon
t  = 0:h:Tf-h

Ms = 1.7 # Maximum allowed magnitude of sensitivity function
Mt = 1.7 # Maximum allowed magnitude of complimentary sensitivity function

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
    J0     = mean(abs,y) # ~ Integrated absolute error IAE
    J = J0
    J -= 0.1*sum(x->log(0.11*max(x,1e-5)),Ms-Sw) # Log-barrier representation of constraint
    J -= 0.1*sum(x->log(0.11*max(x,1e-5)),Mt-Tw) # Log-barrier representation of constraint
    J
end
costfun(p0)

# Evaluates and plots candidate controller
function evalsol(p::Vector)
    C     = K(p)
    S     = 1/(1+P*C) # Sensitivity fun
    PS    = ss(P*S)   # TF from load disturbance to output
    T     = tf(1.) - S# Comp. Sensitivity fun
    s     = Simulator(PS, (t,x) -> [1]) # Sim. unit step load disturbance
    x0    = zeros(PS.nx)
    ty    = eltype(p) # So that all inputs to solve have same numerical type (ForwardDiff.Dual)
    tspan = (ty(0.),ty(Tf))
    sol   = solve(s, ty.(x0), tspan, Tsit5())
    y     = PS.C*sol(t) # y = C*x
    Sw    = vec(bode(S,Ω)[1]) # Freq. domain constraints
    Tw    = vec(bode(T,Ω)[1]) # Freq. domain constraints

    plot(t,y', layout=2, show=false)
    plot!(Ω, [Sw Tw] , lab=["Sw" "Tw"], subplot=2, xscale=:log10, yscale=:log10, show=false)
    plot!([Ω[1],Ω[end]], [Ms,Ms], c = :black, l=:dash, subplot=2, show=false)
    plot!([Ω[1],Ω[end]], [Mt,Mt], c = :purple, l=:dash, subplot=2)
    gui()
end

evalsol(res::Optim.MultivariateOptimizationResults) = evalsol(Optim.minimizer(res))
evalsol(res::BlackBoxOptim.OptimizationResults) = evalsol(best_candidate(res))


# Optimize costfun using Optim with automatic differentiation
using Optim
od = OnceDifferentiable(costfun, p; autodiff = :forward)
@time res3 = optimize(od, p, BFGS(), Optim.Options(show_trace=true))
evalsol(res3)

# Try also second order method
td = TwiceDifferentiable(costfun, p; autodiff = :forward)
@time res4 = optimize(td, p, Newton(), Optim.Options(show_trace=true))

# Optimize costfun using derivative-free method
using BlackBoxOptim
res1 = compare_optimizers(costfun; SearchRange = (0.,2.), NumDimensions = length(p0), MaxTime = 20.0)
res2 = bboptimize(costfun, NumDimensions = length(p0), MaxTime = 20.0)
evalsol(res2)
p = best_candidate(res2)
