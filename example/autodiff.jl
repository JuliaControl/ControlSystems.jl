using ControlSystems, OrdinaryDiffEq, NLopt, BlackBoxOptim, ForwardDiff, Plots
p0          = [0.2,0.8,1] # Initial guess
K(kp,ki,kd) = pid(kp=kp, ki=ki, kd=kd)
K(p)        = K(p...)

# Define process model
ζ  = 0.1
ω  = 1.; ω² = ω^2
P  = tf(ω²,[1, 2ζ*ω, ω²])*tf(1,[1,1])

Ω  = exp10.(LinRange(-1,2,150))  # Frequency vector to eval constraints
Ts  = 0.1 # Sample time for time-domain evaluation
tfinal = 60.  # Time horizon
t  = 0:Ts:tfinal-Ts

Ms = 1.4 # Maximum allowed magnitude of sensitivity function
Mt = 1.4 # Maximum allowed magnitude of complimentary sensitivity function

p  = copy(p0)

function timedomain(p)
    C     = K(p[1], p[2], p[3])
    S     = 1/(1+P*C) # Sensitivity fun
    PS    = ss(P*S)   # TF from load disturbance to output
    s     = Simulator(PS, (t,x) -> [1]) # Sim. unit step load disturbance
    ty    = eltype(p) # So that all inputs to solve have same numerical type (ForwardDiff.Dual)
    x0    = zeros(PS.nx) .|> ty
    tspan = (ty(0.),ty(tfinal))
    sol   = solve(s, x0, tspan, Tsit5())
    y     = PS.C*sol(t) # y = C*x
    C,y
end

function freqdomain(p)
    C     = K(p[1], p[2], p[3])
    S     = 1/(1+P*C) # Sensitivity fun
    T     = tf(1.) - S# Comp. Sensitivity fun
    Sw    = vec(bode(S,Ω)[1]) # Freq. domain constraints
    Tw    = vec(bode(T,Ω)[1]) # Freq. domain constraints
    Sw,Tw
end

# Evaluates and plots candidate controller
function evalsol(p::Vector)
    C,y = timedomain(p)
    Sw,Tw = freqdomain(p)
    plot(t,y', layout=2, show=false)
    plot!(Ω, [Sw Tw] , lab=["Sw" "Tw"], subplot=2, xscale=:log10, yscale=:log10, show=false)
    plot!([Ω[1],Ω[end]], [Ms,Ms], c = :black, l=:dash, subplot=2, show=false, lab="Ms")
    plot!([Ω[1],Ω[end]], [Mt,Mt], c = :purple, l=:dash, subplot=2, lab="Mt")
    gui()
    false
end
evalsol(res::BlackBoxOptim.OptimizationResults) = evalsol(best_candidate(res))


function constraintfun(p)
    Sw,Tw = freqdomain(p)
    [maximum(Sw)-Ms; maximum(Tw)-Mt]
end


function costfun(p)
    C,y = timedomain(p)
    mean(abs,y) # ~ Integrated absolute error IAE
end


function runopt(p, costfun, constraintfun;
    f_tol = 1e-5,
    x_tol = 1e-3,
    c_tol = 1e-8,
    f_cfg = ForwardDiff.GradientConfig(costfun, p),
    g_cfg = ForwardDiff.JacobianConfig(constraintfun, p),
    lb = zeros(length(p)))

    c1 = constraintfun(p)
    np = length(p)
    nc = length(c1)

    function f(p::Vector, grad::Vector)
        if length(grad) > 0
            grad .= ForwardDiff.gradient(costfun,p,f_cfg)
        end
        costfun(p)
    end

    function c(result, p::Vector, grad)
        if length(grad) > 0
            grad .= ForwardDiff.jacobian(constraintfun,p,g_cfg)'
        end
        result .= constraintfun(p)
    end

    opt = Opt(:LD_SLSQP, np)
    lower_bounds!(opt, lb)
    xtol_rel!(opt, x_tol)
    ftol_rel!(opt, f_tol)

    min_objective!(opt, f)
    inequality_constraint!(opt, c, c_tol*ones(nc))
    minf,minx,ret = NLopt.optimize(opt, p)
end

f_cfg = ForwardDiff.GradientConfig(costfun, p)
g_cfg = ForwardDiff.JacobianConfig(constraintfun, p)
@time minf,minx,ret = runopt(1p0, costfun, constraintfun, x_tol=1e-6, c_tol=1e-12, f_cfg=f_cfg, g_cfg=g_cfg)
evalsol(minx)

# # Optimize costfun using derivative-free method
# res1 = compare_optimizers(costfun; SearchRange = (0.,2.), NumDimensions = length(p0), MaxTime = 20.0)
# res2 = bboptimize(costfun, NumDimensions = length(p0), MaxTime = 20.0)
# evalsol(res2)
# p = best_candidate(res2)
