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

This example will demonstrate simple usage of AD using ForwardDiff.jl for optimization-based autotuning of a PID controller.

The system we will control is a double-mass system, in which two masses (or inertias) are connected by a flexible transmission. 

We start by defining the system model and an initial guess for the PID controller parameters
```@example autodiff
using ControlSystemsBase, ForwardDiff, Plots

P = DemoSystems.double_mass_model()

bodeplot(P, title="Bode plot of Double-mass system \$P(s)\$")
```

```@example autodiff
Ω = exp10.(-2:0.04:3)
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

To make the automatic gradient computation through the matrix exponential used in the function [`c2d`](@ref) work, we load the package `ChainRules` that contains a rule for `exp!`, and `ForwardDiffChainRules` that makes ForwardDiff understand the rules in `ChainRules`. Lastly, we need to tell ForwardDiff to use the `exp!` rule for the matrix exponential, which we do by defining an appropriate `@ForwardDiff_frule` for `exp!` that uses the rule in `ChainRules`. All other functions we used work out of the box with ForwardDiff.

```@example autodiff
using Optimization, Statistics, LinearAlgebra
using OptimizationGCMAES
using ChainRules, ForwardDiffChainRules
@ForwardDiff_frule LinearAlgebra.exp!(x1::AbstractMatrix{<:ForwardDiff.Dual})

function plot_optimized(P, params, res)
    fig = plot(layout=(1,3), size=(1200,400), bottommargin=2Plots.mm)
    for (i,params) = enumerate((params, res))
        ls = (i == 1 ? :dash : :solid)
        lab = (i==1 ? "Initial" : "Optimized")
        C, G, r1, r2 = systems(P, params)
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

function systems(P, params)
    kp,ki,kd,Tf = params # We optimize parameters in 
    C    = pid(kp, ki, kd; form=:parallel, Tf, state_space=true)
    G    = extended_gangoffour(P, C) # [S PS; CS T]
    Gd   = c2d(G, 0.1)   # Discretize the system
    res1 = step(Gd[1,1], 0:0.1:15) # Simulate S
    res2 = step(Gd[1,2], 0:0.1:15) # Simulate PS
    C, G, res1, res2
end

σ(x) = 1/(1 + exp(-x)) # Sigmoid function used to obtain a smooth constraint on the peak of the sensitivity function

@views function cost(P, params::AbstractVector{T}) where T
    C, G, res1, res2 = systems(P, params)
    R,_ = bode(G, Ω, unwrap=false)
    S = sum(σ.(100 .* (R[1, 1, :] .- Msc))) # max sensitivity
    CS = sum(Ω .* R[2, 1, :])               # max noise sensitivity
    perf = mean(abs2, res1.y .* res1.t) + mean(abs2, res2.y .* res2.t)
    return perf + Sweight*S + CSweight*CS # Blend all objectives together
end

Msc = 1.3        # Constraint on Ms
Sweight  = 10    # Sensitivity violation penalty
CSweight = 0.001 # Noise amplification penalty

params  = [kp, ki, kd, 0.01] # Initial guess for parameters
using Optimization
using OptimizationGCMAES

fopt = OptimizationFunction((x, _)->cost(P, x), Optimization.AutoForwardDiff())
prob = OptimizationProblem(fopt, params, lb=zeros(length(params)), ub = 10ones(length(params)))
solver = GCMAESOpt()
res = solve(prob, solver; maxiters=1000); res.objective
plot_optimized(P, params, res.u)
```

The optimized controller achieves more or less the same low peak in the sensitivity function, but does this while *both* making the step responses significantly faster *and* using much less controller gain for large frequencies (the orange sensitivity function), an altogether better tuning. The only potentially negative effect of this tuning is that the overshoot in response to a reference step increased slightly, indicated also by the slightly higher peak in the complimentary sensitivity function (green). However, the response to reference steps can (and most often should) be additionally shaped by reference pre-filtering (sometimes referred to as "feedforward" or "reference shaping"), by introducing an additional filter appearing in the feedforward path only, thus allowing elimination of the overshoot without affecting the closed-loop properties.