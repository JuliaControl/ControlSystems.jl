using ControlSystems, OrdinaryDiffEq, Optim, BlackBoxOptim
p0          = Float64[0.2,0.8,2] # Initial guess
K(kp,ki,kd) = pid(kp=kp, ki=ki, kd=kd)
K(p)        = K(p...)

# Define process model
ζ  = 0.1
ω  = 1.; ω² = ω^2
P  = tf(ω²,[1, 2ζ*ω, ω²])*tf(1,[1,1])

Ω  = logspace(-1,2,150)  # Frequency vector to eval constraints
h  = 0.1 # Sample time for time-domain evaluation
Tf = 60.  # Time horizon
t  = 0:h:Tf-h

Ms = 2. # Maximum allowed magnitude of sensitivity function
Mt = 2. # Maximum allowed magnitude of complimentary sensitivity function

p  = copy(p0)

function logbarrier(x, γ=1)
    y = similar(x)
    for i in eachindex(x)
        y[i] = x[i] > 0 ? -log(x[i])/γ : 100000
    end
    y
end

function softmax(x,γ=1)
    mx = maximum(x)
    mx + log(sum(exp.((x.-mx).*γ)))/γ
end

function simulate(p)
    C     = K(p[1], p[2], p[3])
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
    C,y,Sw,Tw
end

# Evaluates and plots candidate controller
function evalsol(p::Vector)
    C,y,Sw,Tw = simulate(p)
    plot(t,y', layout=2, show=false)
    plot!(Ω, [Sw Tw] , lab=["Sw" "Tw"], subplot=2, xscale=:log10, yscale=:log10, show=false)
    plot!([Ω[1],Ω[end]], [Ms,Ms], c = :black, l=:dash, subplot=2, show=false, lab="Ms")
    plot!([Ω[1],Ω[end]], [Mt,Mt], c = :purple, l=:dash, subplot=2, lab="Mt")
    gui()
    false
end
evalsol(res::Optim.MultivariateOptimizationResults) = evalsol(Optim.minimizer(res))
evalsol(state::Optim.OptimizationState) = evalsol(state.metadata["x"])
evalsol(res::BlackBoxOptim.OptimizationResults) = evalsol(best_candidate(res))

function lagrangean(p, λ, s)
    C,y,Sw,Tw = simulate(p)
    L = costfun(y)
    g = constraintfun(Sw,Tw)
    L += λ'*(g-s)
    L += logbarrier([p; s])
end

function constraintfun(Sw,Tw)
    [Ms-softmax(Sw,2000),
    Mt-softmax(Tw,2000)]
end

costfun(y) =   mean(abs,y) # ~ Integrated absolute error IAE

function isfeasible(p)
    C,y,Sw,Tw = simulate(p)
    all(constraintfun(Sw,Tw) .> 0)
end

function opt(p0)
    np = length(p0)
    nc = 2
    C,y,Sw,Tw = simulate(p0)
    pinds = 1:np
    λinds = np+1:np+nc
    sinds = np+nc+1:np+2nc
    if !isfeasible(p0)
        p0 = find_feasible()
        C,y,Sw,Tw = simulate(p)
    end
    s    = constraintfun(Sw,Tw)
    λ    = zeros(size(s))
    p    = [p0; λ; s]
    f(p) = lagrangean(p[pinds], p[λinds], p[sinds])
    cfg  = ForwardDiff.GradientConfig(f, p)
    grad(p) = ForwardDiff.gradient(f,p,cfg)
    α = 0.1

    for i = 1:100
        p .-= α*grad(p)
        @show f(p)
    end
    p[pinds], p[λinds], p[sinds]
end

opt(p0)

# Optimize costfun using Optim with automatic differentiation

# _,_,Sw,Tw = simulate(p0)
# pe0       = [p0; max(1e-5,Ms-softmax(Sw,2000)); 0; max(1e-5,Mt-softmax(Tw,2000)); 0]
# opts = Optim.Options(show_trace=true) # Uncomment to turn off live plotting
# lower = [0., 0, 0, 0, -Inf, 0, -Inf]
# upper = [Inf, Inf, Inf, Inf, Inf, Inf, Inf]
# barrier = 100.
# od    = OnceDifferentiable(costfun, pe0; autodiff = :forward)
# res3  = optimize(od, pe0, lower, upper, Fminbox{GradientDescent}(); optimizer_o = Optim.Options(show_trace=false, iterations = 5), show_trace=true, callback=callback)
# evalsol(res3)


# opts      = Optim.Options(show_trace=true, extended_trace=true, callback=callback)
# od    = OnceDifferentiable(costfun, p0; autodiff = :forward)
# @time res3 = optimize(od, p0, GradientDescent(), opts)

# # Try also second order method
# td = TwiceDifferentiable(costfun, p; autodiff = :forward)
# @time res4 = optimize(td, p, NewtonTrustRegion(), opts)

# # Optimize costfun using derivative-free method
# res1 = compare_optimizers(costfun; SearchRange = (0.,2.), NumDimensions = length(p0), MaxTime = 20.0)
# res2 = bboptimize(costfun, NumDimensions = length(p0), MaxTime = 20.0)
# evalsol(res2)
# p = best_candidate(res2)
