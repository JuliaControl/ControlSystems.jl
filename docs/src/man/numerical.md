# Performance considerations

## Numerical accuracy
Transfer functions, and indeed polynomials in general, are infamous for having poor numerical properties. Consider the simple polynomial $ax^n - 1$ which, due to rounding of the polynomial coefficients, is represented as $(a+\epsilon)x^n - 1$ where $\epsilon$ is on the order of `eps(a)`. The roots of this polynomial have a much larger $\epsilon$, due to the n:th root in the expression $\dfrac{1}{\sqrt[n]{(a + \epsilon)}}$. For this reason, it's ill-advised to use high-order transfer functions. Orders as low as 6 may already be considered high. When a transfer function is converted to a state-space representation using `ss(G)`, balancing is automatically performed in an attempt at making the numerical properties of the model better.



## Frequency-response calculation
For small systems (small number of inputs, outputs and states), evaluating the frequency-response of a transfer function is reasonably accurate and very fast.

```julia
G = tf(1, [1, 1])
w = exp10.(LinRange(-2, 2, 200));
@btime freqresp($G, $w);
# 4.351 μs (2 allocations: 3.31 KiB)
```
Evaluating the frequency-response for the equivalent state-space system incurs some additional allocations due to a Hessenberg matrix factorization

```julia
sys = ss(G);
@btime freqresp($sys, $w);
# 20.820 μs (16 allocations: 37.20 KiB)
```

For larger systems, the state-space calculations are considerably more accurate, provided that the realization is well balanced.

For optimal performance, one may preallocate the return array
```julia
ny,nu = size(G)
R = zeros(ComplexF64, ny, nu, length(w));

@btime freqresp!($R, $G, $w);
# 4.214 μs (1 allocation: 64 bytes)
```

Other functions that accept preallocated workspaces are
- [`bodemag!`](@ref)
- [`freqresp!`](@ref)
- [`lsim!`](@ref)

an example using [`bodemag!`](@ref) follows:
```julia
using ControlSystemsBase
G = tf(ssrand(2,2,5))
w = exp10.(LinRange(-2, 2, 20000))
@btime bode($G, $w);
# 55.120 ms (517957 allocations: 24.42 MiB)
@btime bode($G, $w, unwrap=false); # phase unwrapping is slow
# 3.624 ms (7 allocations: 2.44 MiB)
ws = ControlSystemsBase.BodemagWorkspace(G, w)
@btime bodemag!($ws, $G, $w);
# 2.991 ms (1 allocation: 64 bytes)
```



## Time-domain simulation

### Time scale
When simulating a dynamical system in continuous time, a differential-equation integrator is used. These integrators are sensitive to the scaling of the equations, and may perform poorly for stiff problems or problems with a poorly chosen time scale. In, e.g., electronics, it's common to simulate systems where the dominant dynamics have time constants on the order of microseconds. To simulate such systems accurately, it's often a good idea to model the system in microseconds rather than in seconds. The function [`time_scale`](@ref) can be used to automatically change the time scale of a system.

### Transfer functions
Transfer functions are automatically converted to state-space form before time-domain simulation. If you want control over the exact internal representation used, consider modeling the system as a state-space system already from start. 

### Discrete-time simulation
Linear systems with zero-order-hold inputs can be exactly simulated in discrete time. You may specify ZoH-discretization in the call to [`lsim`](@ref) using `method=:zoh` or manually perform the discretization using [`c2d`](@ref). Discrete-time simulation is often *much* faster than continuous-time integration.

For discrete-time systems, the function [`lsim!`](@ref) accepts a pre-allocated workspace objects that can be used to avoid allocations for repeated simulations.

