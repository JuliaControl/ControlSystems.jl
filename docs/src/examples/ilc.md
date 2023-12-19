# Iterative-Learning Control
In this example, we will design an [Iterative-Learning Control (ILC)](https://en.wikipedia.org/wiki/Iterative_learning_control) iteration scheme. ILC can be though of as a simple reinforcement-learning strategy that is suitable in situations where a *repetitive task* is to be performed multiple times, and disturbances acting on the system are also repetitive and predictable but unknown. Multiple versions of ILC exists, in this tutorial we will consider a heuristic scheme as well as a model-based scheme. 

## Algorithm

The ILC iteration scheme typically looks something like this (many variants exists), at ILC iteration $k$:
```math
\begin{aligned}
y_k(t) &= G(q) \big(r(t) + a_k(t) \big) \\
e_k(t) &= r(t) - y_k(t) \\
a_k(t) &= Q(q) \big( a_{k-1}(t) + L(q) e_{k-1}(t) \big)
\end{aligned}
```
where $q$ is the time-shift operator, $G(q)$ is the transfer function from the reference $r$ to the output $y$, i.e, typically a closed-loop transfer function, $e_k$ is the control error and $a_k$ is the ILC adjustment signal, an additive correction to the reference that is learned throughout the ILC iterations in order to minimize the control error. $Q(q)$ and $L(q)$ are stable filters that control the learning dynamics. Interestingly, these filters does not have to be causal since they operate on the signals $e$ and $a$ *between* ILC iterations, when the whole signals are available at once for acausal filtering. 

In simulation (the rollout $y_k = G(q) (r + a_k)$ is simulated), this scheme is nothing other than an open-loop optimal-control strategy, while if $y_k = G(q) (r + a_k)$ amounts to performing an actual experiment on a process, ILC turns into episode-based reinforcement learning or adaptive control.

The system to control in this example is a double-mass system with a spring and damper in between. This system is a common model of a servo system where one mass represents the motor and the other represents the load. The spring and damper represents a flexible transmission between them. We will create two instances of the system model. ``G`` represents the nominal model, whereas ``G_{act}`` represents the actual (unknown) dynamics. This simulates a model-based approach where there is a slight error in the model. The error will lie in the mass of the load, simulating, e.g., that the motor is driving a heavier load than specified. 

## System model and controller

```@example ilc
using ControlSystemsBase, Plots

function double_mass_model(; 
                Jm = 1,   # motor inertia
                Jl = 1,   # load inertia
                k  = 100, # stiffness
                c0 = 1,   # motor damping
                c1 = 1,   # transmission damping
                c2 = 1,   # load damping
)

    A = [
        0.0 1 0 0
        -k/Jm -(c1 + c0)/Jm k/Jm c1/Jm
        0 0 0 1
        k/Jl c1/Jl -k/Jl -(c1 + c2)/Jl
    ]
    B = [0, 1/Jm, 0, 0]
    C = [1 0 0 0]
    ss(A, B, C, 0)
end

G    = double_mass_model(Jl = 1)
Gact = double_mass_model(Jl = 1.5) # 50% more load than modeled

bodeplot([G, Gact], lab=["G model" "G actual"], plotphase=false)
```
We will design a PID controller with a filter for the system, the controller is poorly tuned and not very good at tracking fast reference steps, in practice, one would likely design a feedforward controller as well to improve upon this, but for now we'll stick with the simple feedback controller.

```@example ilc
C  = pid(10, 1, 1, form = :series) * tf(1, [0.02, 1])
Ts = 0.02 # Sample time
Gc = c2d(feedback(G*C), Ts)       |> tf
Gcact = c2d(feedback(Gact*C), Ts) |> tf
plot(step(Gc, 10), title="Closed-loop step response", lab="model")
plot!(step(Gcact, 10), lab="actual")
```

## Reference trajectory

Next up we design a reference trajectory and simulate the actual closed-loop dynamics.
```@example ilc
T = 3pi    # Duration
t = 0:Ts:T # Time vector
function funnysin(x)
    x = sin(x)
    s,a = sign(x), abs(x)
    s*((a + 0.01)^0.2 - 0.01^0.2)
end
r = funnysin.(t)' |> Array # Reference signal

res = lsim(Gcact, r, t)
plot(res, plotu=true, layout=1, sp=1, title="Closed-loop simulation with actual dynamics", lab=["y" "r"])
```
Performance is poor.. Enter ILC!

## Non-causal filtering

For ILC to work well, we define two helper functions. One that applies a zero-phase filter by filtering both forwards and backwards (`filtfilt`). This is possible since ILC operates on signals offline, between iterations in the ILC scheme. We also define a special `lsim` that handles non-causal systems to allow "lookahead" into the future. This typically improves the performance of ILC by quite a lot, and is once again possible since ILC operates on prerecorded signals. 

```@example ilc
function lsim_zerophase(G, u, args...; kwargs...)
    res = lsim(G, u[:, end:-1:1], args...; kwargs...)
    lsim(G, res.y[:, end:-1:1], args...; kwargs...).y
end

function lsim_noncausal(L::LTISystem{<:Discrete}, u, args...; kwargs...)
    np = length(denpoly(L)[])
    nz = length(numpoly(L)[])
    zeroexcess = nz-np
    if zeroexcess <= 0
        return lsim(L, u, args...; kwargs...)
    end
    integrators = tf(1, [1, 0], L.Ts)^zeroexcess
    res = lsim(L*integrators, u, args...; kwargs...)
    res.y[1:end-zeroexcess] .= res.y[1+zeroexcess:end]
    res.y
end
nothing # hide
```

## Choosing filters
The next step is to define the ILC filters ``Q(x)`` and ``L(z)``.

The filter $L(q)$ acts as a frequency-dependent step size. To make the procedure take smaller steps, simply scale $L$ by a constant < 1. Scaling down $L$ makes the learning process slower but more robust. A heuristic choice of $L$ is some form of scaled lookahead, such as $0.5z^l$ where $l \geq 0$ is the number of samples lookahead. A model-based approach may use some form of inverse of the system model, which is what we will use here. [^nonlinear]

[^nonlinear]: Inverse models can be formed also for some nonlinear systems. [ModelingToolkit.jl](https://mtk.sciml.ai/dev/) is particularily well suited for inverting models due to its acausal nature.

The filter $Q(q)$ acts to make the procedure robust w.r.t. noise and modeling errors. $Q$ has a final say over what frequencies appear in $a$ and it's good to choose $Q$ with low-pass properties. $Q$ will here be applied in zero-phase mode, so the effective transfer function will be $Q(z)Q(z̄)$.

```@example ilc
z = tf("z", Ts)
Q = c2d(tf(1, [0.05, 1]), Ts)
# L = 0.9z^1 # A more conservative and heuristic choice
L = 0.5inv(Gc) # Make the scaling factor smaller to take smaller steps
nothing # hide
```

A theorem due to Norrlöf says that for the ILC iterations to converge, one needs to satisfy
$$| 1 - LG | < |Q^{-1}|$$
which we can verify by looking at the Bode curves of the two sides of the inequality
```@example ilc
bodeplot([inv(Q), (1 - L*Gc)], plotphase=false, lab=["Stability boundary \$Q^{-1}\$" "\$1 - LG\$"])
bodeplot!((1 - L*Gcact), plotphase=false, lab="\$1 - LG\$ actual")
```
Above, we plotted this curve also for the actual dynamics. This is of course not possible in a real scenario where this is unknown, but one could plot it for multiple plausible models and verify that they are all below the boundary. See [Uncertainty modeling using RobustAndOptimalControl.jl](https://juliacontrol.github.io/RobustAndOptimalControl.jl/dev/uncertainty/) for guidance on this. Looking at the stability condition, it becomes obvious how making $Q$ small where the model is uncertain is beneficial for robustness of the ILC scheme.

## ILC iteration

The next step is to implement the ILC scheme and run it:

```@example ilc
function ilc(Gc, Q, L)
    a = zero(r) # ILC adjustment signal starts at 0
    fig1 = plot(t, vec(r), sp=1, layout=(3,1), l=(:black, 3), lab="Ref")
    fig2 = plot(title="Sum of squared error", xlabel="Iteration", legend=false, titlefontsize=10, framestyle=:zerolines, ylims=(0, 7.1))
    for iter = 1:5
        ra = r .+ a
        res = lsim(Gc, ra, t) # Simulate system, replaced by an actual experiment if running on real process
        y = res.y             # System response
        e = r .- y            # Error
        Le = lsim_noncausal(L, e, t)
        a  = lsim_zerophase(Q, a + Le, t) # Update ILC adjustment

        err = sum(abs2, e)
        plot!(fig1, res, plotu=true, sp=[1 2], title=["Output \$y(t)\$" "Adjusted reference \$r + a\$"], lab="Iter $iter", c=iter)
        plot!(fig1, e[:], sp=3, title="Tracking error \$e(t)\$", lab="err: $(round(err, digits=2))", c=iter)
        scatter!(fig2, [iter], [err])
    end
    plot(fig1, fig2, layout=@layout([a{0.7w} b{0.3w}]))
end
ilc(Gc, Q, L)
```
When running on the model, the result looks very good.
We see that the tracking error in the last plot decreases rapidly and is much smaller after only a couple of iterations. We also note that the adjusted reference $r+a$ has effectively been phase-advanced slightly to compensate for the lag in the system dynamics. This is an effect of the acausal filtering due to $L = G_C^{-1}$.


How does it work on the "actual" dynamics?
```@example ilc
ilc(Gcact, Q, L)
```
The result is subtly worse, but considering the rather big model error the result is still quite good. 
## Summary
We have seen how ILC can be used to improve tracking performance in a scenario where a repetitive task is to be executed several times. In simulation like here, ILC can be seen as an optimal-control strategy to come up with a optimal reference trajectory to minimize the control error, while if implemented on a physical process, the scheme amounts to a simple but effective reinforcement-learning or adaptive-control approach. ILC often works very well in practice and has been used in robotics and machining among other areas. 

ILC does not work very well if stochastic disturbances dictate the control performance or a task is to be performed only a small number of times. In, e.g., machining applications, each ILC iteration may imply performing destructive machining on expensive material with suboptimal result before convergence. This may only be cost effective if the task is to be performed many times after an initial "tuning" by means of ILC.
