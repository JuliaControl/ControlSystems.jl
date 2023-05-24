# Automatic Differentiation
In Julia, it is often possible to automatically compute derivatives, gradients, Jacobians and Hessians of arbitrary Julia functions with precision matching the machine precision, that is, without the numerical inaccuracies incurred by finite-difference approximations.

Two general methods for automatic differentiation are available: forward and reverse mode. Forward mode is algorithmically more favorable for functions with few inputs but many outputs, while reverse mode is more efficient for functions with many parameters but few outputs (like in deep learning). In Julia, forward-mode AD is provided by the package [ForwardDiff.jl](https://github.com/JuliaDiff/ForwardDiff.jl), while reverse-mode AD is provided by several different packages, such as Zygote.jl and ReverseDiff.jl. Forward-mode AD generally has a lower overhead than reverse-mode AD, so for functions of a small number of parameters, say, less than about 10 or 100, forward-mode is usually most efficient. ForwardDiff.jl also has support for differentiating most of the Julia language, making the probability of success higher than for other packages, why we generally recommend trying ForwardDiff.jl first.

## Linearizing nonlinear dynamics
Nonlinear dynamics on the form
```math
\begin{aligned}
\dot x &= f(x, u) \\
y &= g(x, u)
\end{aligned}
```

is easily linearized in the point ``x_0, u_0`` using ForwardDiff.jl:
```@example autodiff
using ControlSystemsBase, ForwardDiff

"An example of nonlinear dynamics"
function f(x, u)
    x1, x2 = x
    u1, u2 = u
    [x2; u1*x1 + u2*x2]
end

x0 = [1.0, 0.0] # Operating point to linearize around
u0 = [0.0, 1.0]

A = ForwardDiff.jacobian(x -> f(x, u0), x0)
B = ForwardDiff.jacobian(u -> f(x0, u), u0)

"An example of a nonlinear output (measurement) function"
function g(x, u)
    y = [x[1] + 0.1x[1]*u[2]; x[2]]
end

C = ForwardDiff.jacobian(x -> g(x, u0), x0)
D = ForwardDiff.jacobian(u -> g(x0, u), u0)

linear_sys = ss(A, B, C, D)
```
The example above linearizes `f` in the point ``x_0, u_0`` to obtain the linear statespace matrices ``A`` and ``B``, and linearizes `g` to obtain the linear output matrices ``C`` and ``D``.

## Optimization-based tuning

### PID controller

This example will demonstrate simple usage of AD using ForwardDiff.jl for optimization-based auto tuning of a PID controller.

The system we will control is a double-mass system, in which two masses (or inertias) are connected by a flexible transmission. 

We start by defining the system model and an initial guess for the PID controller parameters
```@example autodiff
using ControlSystemsBase, ForwardDiff, Plots

P = DemoSystems.double_mass_model()

bodeplot(P, title="Bode plot of Double-mass system \$P(s)\$")
```

```@example autodiff
Ω = exp10.(-2:0.03:3)
kp,ki,kd,Tf =  1, 0.1, 0.1, 0.01 # controller parameters

C  = pid(kp, ki, kd; Tf, form=:parallel, state_space=true) # Construct a PID controller with filter
G  = feedback(P*C) # Closed-loop system
S  = 1/(1 + P*C)   # Sensitivity function
Gd = c2d(G, 0.1)   # Discretize the system
res = step(Gd,15)  # Step-response

mag = bodev(S, Ω)[1]
plot(res, title="Time response", layout = (1,3), legend=:bottomright)
plot!(Ω, mag, title="Sensitivity function", xscale=:log10, yscale=:log10, subplot=2, legend=:bottomright, ylims=(3e-2, Inf))
Ms, _ = hinfnorm(S)
hline!([Ms], l=(:black, :dash), subplot=2, lab="\$M_S = \$ $(round(Ms, digits=3))", sp=2)
nyquistplot!(P*C, Ω, sp=3, ylims=(-2.1,1.1), xlims=(-2.1,1.2), size=(1200,400))
```
The initial controller ``C`` achieves a maximum peak of the sensitivity function of ``M_S = 1.3`` which implies a rather robust tuning, but the step response is sluggish. We will now try to optimize the controller parameters to achieve a better performance.

We start by defining a helper function `plot_optimized` that will evaluate the performance of the tuned controller. We then define a function `systems` that constructs the gang-of-four transfer functions ([`extended_gangoffour`](@ref)) and performs time-domain simulations of the transfer functions ``S(s)`` and ``P(s)S(s)``, i.e., the transfer functions from reference ``r`` to control error ``e``, and the transfer function from an input load disturbance ``d`` to the control error ``e``. By optimizing these step responses with respect to the PID parameters, we will get a controller that achieves good performance. To promote robustness of the closed loop as well as to limit the amplification of measurement noise in the control signal, we penalize the peak of the sensitivity function ``S`` as well as the (approximate) frequency-weighted ``H_2`` norm of the transfer function ``CS(s)``.


The cost function `cost` encodes the constraint on the peak of the sensitivity function as a penalty function (this could be enforced explicitly by using a constrained optimizer) and weighs the different objective terms together using user-defined weights `Sweight` and `CSweight`. Finally, we use [Optimization.jl](https://github.com/SciML/Optimization.jl) to optimize the cost function and tell it to use ForwardDiff.jl to compute the gradient of the cost function. The optimizer we use in this example, `GCMAESOpt`, is the "Gradient-based Covariance Matrix Adaptation Evolutionary Strategy", which can be thought of as a blend between a derivative-free global optimizer and a gradient-based local optimizer.

To make the automatic gradient computation through the matrix exponential used in the function [`c2d`](@ref)[^zoh] work, we load the package `ChainRules` that contains a rule for `exp!`, and `ForwardDiffChainRules` that makes ForwardDiff understand the rules in `ChainRules`. Lastly, we need to tell ForwardDiff to use the `exp!` rule for the matrix exponential, which we do by defining an appropriate `@ForwardDiff_frule` for `exp!` that uses the rule in `ChainRules`. All other functions we used work out of the box with ForwardDiff.

[^zoh]: Only applies for the methods that rely on matrix exponential, such as zero and first-order hold, if `c2d(P, Ts, :tustin)` is used instead, `ChainRules` and `ForwardDiffChainRules` are not required.

```@example autodiff
using Optimization, Statistics, LinearAlgebra
using OptimizationGCMAES, ChainRules, ForwardDiffChainRules
# @ForwardDiff_frule LinearAlgebra.exp!(x1::AbstractMatrix{<:ForwardDiff.Dual})

function plot_optimized(P, params, res)
    fig = plot(layout=(1,3), size=(1200,400), bottommargin=2Plots.mm)
    for (i,params) = enumerate((params, res))
        ls = (i == 1 ? :dash : :solid)
        lab = (i==1 ? "Initial" : "Optimized")
        C, G, r1, r2 = systems(params, P)
        mag = reshape(bode(G, Ω)[1], 4, :)'[:, [1, 2, 4]]
        plot!([r1, r2]; title="Time response", subplot=1,
            lab = lab .* [" \$r → e\$" " \$d → e\$"], legend=:bottomright, ls,
            fillalpha=0.05, linealpha=0.8, seriestype=:path, c=[1 3])
        plot!(Ω, mag; title="Sensitivity functions \$S(s), CS(s), T(s)\$",
            xscale=:log10, yscale=:log10, subplot=2, lab, ls,
            legend=:bottomright, fillalpha=0.05, linealpha=0.8, c=[1 2 3], linewidth=i)
        nyquistplot!(P*C, Ω; Ms_circles=Msc, sp=3, ylims=(-2.1,1.1), xlims=(-2.1,1.2), lab, seriescolor=i, ls)
    end
    hline!([Msc], l=:dashdot, c=1, subplot=2, lab="Constraint", ylims=(9e-2, Inf))
    fig
end

function systems(params, P)
    kp,ki,kd,Tf = params # We optimize parameters in 
    C    = pid(kp, ki, kd; form=:parallel, Tf, state_space=true)
    G    = extended_gangoffour(P, C) # [S PS; CS T]
    Gd   = c2d(G, 0.1)   # Discretize the system
    res1 = step(Gd[1,1], 0:0.1:15) # Simulate S
    res2 = step(Gd[1,2], 0:0.1:15) # Simulate PS
    C, G, res1, res2
end

σ(x) = 1/(1 + exp(-x)) # Sigmoid function used to obtain a smooth constraint on the peak of the sensitivity function

@views function cost(params::AbstractVector{T}, P) where T
    C, G, res1, res2 = systems(params, P)
    R,_ = bode(G, Ω, unwrap=false)
    S = sum(σ.(100 .* (R[1, 1, :] .- Msc))) # max sensitivity
    CS = sum(Ω .* R[2, 1, :])               # frequency-weighted noise sensitivity
    perf = mean(abs2, res1.y .* res1.t') + mean(abs2, res2.y .* res2.t')
    return perf + Sweight*S + CSweight*CS # Blend all objectives together
end

Msc = 1.3        # Constraint on Ms
Sweight  = 10    # Sensitivity violation penalty
CSweight = 0.001 # Noise amplification penalty

params  = [kp, ki, kd, 0.01] # Initial guess for parameters
using Optimization
using OptimizationGCMAES

fopt = OptimizationFunction(cost, Optimization.AutoForwardDiff())
prob = OptimizationProblem(fopt, params, P, lb=zeros(length(params)), ub = 10ones(length(params)))
solver = GCMAESOpt()
res = solve(prob, solver; maxiters=1000)
plot_optimized(P, params, res.u)
```

The optimized controller achieves more or less the same low peak in the sensitivity function, but does this while *both* making the step responses significantly faster *and* using much less controller gain for large frequencies (the orange sensitivity function), an altogether better tuning. The only potentially negative effect of this tuning is that the overshoot in response to a reference step increased slightly, indicated also by the slightly higher peak in the complimentary sensitivity function (green). However, the response to reference steps can (and most often should) be additionally shaped by reference pre-filtering (sometimes referred to as "feedforward" or "reference shaping"), by introducing an additional filter appearing in the feedforward path only, thus allowing elimination of the overshoot without affecting the closed-loop properties.

### LQG controller
We could attempt a similar automatic tuning of an LQG controller. This time, we choose to optimize the weight matrices of the LQR problem and the state covariance matrix of the noise. The synthesis of an LQR controller involves the solution of a Ricatti equation, which in turn involves performing a Schur decomposition. These steps hard hard to differentiate through in a conventional way, but we can make use of implicit differentiation using the implicit function theorem. To do so, we load the package `ImplicitDifferentiation`, and define the conditions that hold at the solution of the Ricatti equaiton:
```math
A^TX + XA - XBR^{-1}B^T X + Q = 0
```

When `ImplicitDifferentiation` is loaded, differentiable versions of [`lqr`](@ref) and [`kalman`](@ref) that make use of the "implicit function" are automatically loaded.

```@example autodiff
using ImplicitDifferentiation, ComponentArrays # Both these packages are required to load the implicit differentiation rules
```

Since this is a SISO system, we do not need to tune the control-input matrix or the measurement covariance matrix, any non-unit weight assigned to those can be associated with the state matrices instead. Since these matrices are supposed to be positive semi-definite, we optimize Cholesky factors rather than the full matrices. We will also reformulate the problem slightly and use proper constraints to limit the sensitivity peak.
```@example autodiff
using Ipopt, OptimizationMOI
MOI = OptimizationMOI.MOI
function triangular(x)
    m = length(x)
    n = round(Int, sqrt(2m-1))
    T = zeros(eltype(x), n, n)
    k = 1
    for i = 1:n, j = i:n
        T[i,j] = x[k]
        k += 1
    end
    T
end
invtriangular(T) = [T[i,j] for i = 1:size(T,1) for j = i:size(T,1)]

function systemslqr(params::AbstractVector{T}, P) where T
    n2 = length(params) ÷ 2
    Qchol = triangular(params[1:n2])
    Rchol = triangular(params[n2+1:2n2])
    Q = Qchol'Qchol
    R = Rchol'Rchol
    L = lqr(P, Q, one(T)*I(1))        # The last matrix must have the correct type
    K = kalman(P, R, one(T)*I(1))     # The last matrix must have the correct type
    C = observer_controller(P, L, K)
    G = extended_gangoffour(P, C) # [S PS; CS T]
    C, G
end

function systems(params, P)
    C, G = systemslqr(params, P)
    C, G, sim(G)...
end

function sim(G)
    Gd = c2d(G, 0.1, :zoh)          # Discretize the system
    res1 = step(Gd[1, 1], 0:0.1:15) # Simulate S
    res2 = step(Gd[1, 2], 0:0.1:15) # Simulate PS
    res1, res2
end

@views function cost(params::AbstractVector{T}, P) where T
    CSweight = 0.001 # Noise amplification penalty
    C, G = systemslqr(params, P)
    res1, res2 = sim(G)
    R, _ = bodev(G[2, 1], Ω; unwrap=false)
    CS = sum(R .*= Ω)               # frequency-weighted noise sensitivity
    perf = mean(abs2, res1.y .*= res1.t') + mean(abs2, res2.y .*= res2.t')
    return perf + CSweight * CS # Blend all objectives together
end

@views function constraints(res, params::AbstractVector{T}, P) where T
    C, G = systemslqr(params, P)
    S, _ = bodev(G[1, 1], Ω; unwrap=false)
    res .= maximum(S) # max sensitivity
    nothing
end

Msc = 1.3   # Constraint on Ms

Q0 = diagm([1.0, 1, 1, 1]) # Initial guess LQR state penalty
R0 = diagm([1.0, 1, 1, 1]) # Initial guess Kalman state covariance
params = [invtriangular(cholesky(Q0).U); invtriangular(cholesky(R0).U)]

fopt = OptimizationFunction(cost, Optimization.AutoForwardDiff(); cons=constraints)
prob = OptimizationProblem(fopt, params, P;
    lb    = fill(-10.0, length(params)),
    ub    = fill(10.0, length(params)),
    ucons = fill(Msc, 1),
    lcons = fill(-Inf, 1),
)
solver = Ipopt.Optimizer()
MOI.set(solver, MOI.RawOptimizerAttribute("print_level"), 0)
MOI.set(solver, MOI.RawOptimizerAttribute("acceptable_tol"), 1e-1)
MOI.set(solver, MOI.RawOptimizerAttribute("acceptable_constr_viol_tol"), 1e-2)
MOI.set(solver, MOI.RawOptimizerAttribute("acceptable_iter"), 5)
MOI.set(solver, MOI.RawOptimizerAttribute("hessian_approximation"), "limited-memory")

res = solve(prob, solver)
plot_optimized(P, params, res.u)
```

This controller should perform better than the PID controller, which is known to be incapable of properly damping the resonance in a double-mass system. However, we did not include any integral action in the LQG controller, which has implication for the disturbance response, as indicated by the green step response in the simulation above.


#### Robustness analysis
To check the robustness of the designed LQG controller w.r.t. parametric uncertainty in the plant, we load the package [`MonteCarloMeasurements`](https://github.com/baggepinnen/MonteCarloMeasurements.jl) and recreate the plant model with 20% uncertainty in the spring coefficient.
```@example autodiff
using MonteCarloMeasurements
Pu = DemoSystems.double_mass_model(k = Particles(32, Uniform(80, 120))) # Create a model with uncertainty in spring stiffness k ~ U(80, 120)
unsafe_comparisons(true) # For the Bode plot to work

C,_ = systemslqr(res.u, P)             # Get the controller assuming P without uncertainty
Gu = extended_gangoffour(Pu, C)     # Form the gang-of-four with uncertainty
w = exp10.(LinRange(-1.5, 2, 500))
bodeplot(Gu, w, plotphase=false, ri=false, N=32, ylims=(1e-1, 30), layout=1, sp=1, c=[1 2 4 3], lab=["S" "CS" "PS" "T"])
hline!([Msc], l=:dashdot, c=1, lab="Constraint", ylims=(9e-2, Inf))
```
The uncertainty in the spring stiffness caused an uncertainty in the resonant peak in the sensitivity functions, it's a good thing that we designed a controller that was conservative with a large margin (small ``M_S``) so that all the plausible variations of the plant are expected to behave reasonably well:
```@example autodiff
Gd   = c2d(Gu, 0.05)   # Discretize the system
r1 = step(Gd[1,1], 0:0.05:15) # Simulate S
r2 = step(Gd[1,2], 0:0.05:15) # Simulate PS
plot([r1, r2]; title="Time response",
            lab = [" \$r → e\$" " \$d → e\$"], legend=:bottomright,
            fillalpha=0.05, linealpha=0.8, seriestype=:path, c=[1 3], ri=false, N=32)
```

## Known limitations
The following issues are currently known to exist when using AD through ControlSystems.jl:

### ForwardDiff
[ForwardDiff.jl](https://github.com/JuliaDiff/ForwardDiff.jl) works for a lot of workflows without any intervention required from the user. The following known limitations exist:
- The function `svdvals` does not have a forward rule defined. This means that the functions [`sigma`](@ref) and `opnorm` will not work for MIMO systems with ForwardDiff. SISO, MISO and SIMO systems will, however, work.
- [`hinfnorm`](@ref) requires ImplicitDifferentiation.jl and ComponentArrays.jl to be manually loaded by the user, after which there are implicit differentiation rules defined for [`hinfnorm`](@ref). The implicit rule calls `opnorm`, and is thus affected by the first limitation above for MIMO systems. [`hinfnorm`](@ref) has a reverse rule defined in RobustAndOptimalControl.jl, which is not affected by this limitation.
- [`are`](@ref), [`lqr`](@ref) and [`kalman`](@ref) all require ImplicitDifferentiation.jl and ComponentArrays.jl to be manually loaded by the user, after which there are implicit differentiation rules defined. To invoke the correct method of these functions, it is important that the second matrix (corresponding to input or measurement) has the `Dual` number type, i.e., the `R` matrix in `lqr(P, Q, R)` or `lqr(Continuous, A, B, Q, R)`
- The `schur` factorization is not amenable to differentiation using ForwardDiff. This is the fundamental reason for requireing ImplicitDifferentiation.jl to differentiate through the Ricatti equation solver. `schur` is called in several additional places, including [`balreal`](@ref) and all [`lyap`](@ref) solvers. To make `schur` differentiable, an implicit differentiation rule would be required.

### Reverse-mode AD
- Zygote does not work very well at all, due to
    - Frequent use of mutation for performance
    - Try catch blocks