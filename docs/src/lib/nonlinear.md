# Nonlinear functionality

!!! danger "Experimental"
    The nonlinear interface is currently experimental and at any time subject to breaking changes not respecting semantic versioning. 


ControlSystems.jl can represent nonlinear feedback systems that can be written on the form
```
      ┌─────────┐
 y◄───┤         │◄────u
      │    P    │
Δy┌───┤         │◄───┐Δu
  │   └─────────┘    │
  │                  │
  │      ┌───┐       │
  └─────►│ f ├───────┘
         └───┘
```
i.e., as a linear-fractional transform (LFT) between a linear system ``P`` and a diagonal matrix with scalar non-linear functions $f$. This representation is identical to that used for delay systems, and is exposed to the user in a similar way as well. The main entry point is the function [`nonlinearity`](@ref) which takes a nonlinear function $f$ like so, `nonlinearity(f)`. This creates a primitive system containing only the nonlinearity, but which behaves like a standard `LTISystem` during algebraic operations. We illustrate its usage through a number of examples.

## Examples
### Control-signal saturation
To create a controller that saturates the output at ``\pm 0.7``, we call
```@example nonlinear
using ControlSystems, Plots
using ControlSystemsBase: nonlinearity # This functionality is not exported due to the beta status

C    = pid(1, 0.1, form=:parallel)                  # A standard PI controller
nl   = nonlinearity(x->clamp(x, -0.7, 0.7)) # a saturating nonlinearity
satC = nl*C # Connect the saturation at the output of C
```
we may now use this controller like we would normally do in ControlSystems, e.g.,
```@example nonlinear
P   = tf(1, [1, 1])    # a plant
G   = feedback(P*C)    # closed loop without nonlinearity
Gnl = feedback(P*satC) # closed loop with saturation

Gu   = feedback(C, P)    # closed loop from reference to control signal without nonlinearity
Gunl = feedback(satC, P) # closed loop from reference to control signal with saturation

plot(step([G; Gu], 5), lab = ["Linear y" "Linear u"])
plot!(step([Gnl; Gunl], 5), lab = ["Nonlinear y" "Nonlinear u"])
```

Since the saturating nonlinearity is common, we provide the constructor [`ControlSystemsBase.saturation`](@ref) that automatically forms the equivalent to `nonlinearity(x->clamp(x, -0.7, 0.7))` while at the same time making sure the function has a recognizable name when the system is printed
```@example nonlinear
using ControlSystemsBase: saturation
saturation(0.7)
```

See also [`ControlSystemsBase.ratelimit`](@ref) that saturates the derivative of a signal.

### Non-zero operating point
It's common to linearize nonlinear systems around some operating point. We may make use of the helper constructor [`ControlSystemsBase.offset`](@ref) to create affine functions at the inputs and outputs of the linearized system to, e.g.,
1. Make sure that simulations result are given in the original coordinates rather than in the coordinates of the linearization point.
2. Allow nonlinearities that are added back after the linearization (such as saturations) to operate with their original parameters.

We will demonstrate a composite usage of `offset` and `saturation` below. The system we'll consider is a linearized model of a quadruple-tank process;

The system is linearized around the operating point
```@example nonlinear
xr = [10, 10, 4.9, 4.9] # reference state
ur = [0.263, 0.263]     # control input at the operating point
nothing # hide
```
and is given by
```@example nonlinear
using LinearAlgebra
kc, k1, k2, g = 0.5, 1.6, 1.6, 9.81
A1 = A3 = A2 = A4 = 4.9
a1, a3, a2, a4 = 0.03, 0.03, 0.03, 0.03
h01, h02, h03, h04 = xr
T1, T2 = (A1/a1)sqrt(2*h01/g), (A2/a2)sqrt(2*h02/g)
T3, T4 = (A3/a3)sqrt(2*h03/g), (A4/a4)sqrt(2*h04/g)
c1, c2 = (T1*k1*kc/A1), (T2*k2*kc/A2)
γ1, γ2 = 0.3, 0.3

# Define the process dynamics
A = [-1/T1     0 A3/(A1*T3)          0
        0     -1/T2          0 A4/(A2*T4)
        0         0      -1/T3          0
        0         0          0      -1/T4]
B = [γ1*k1/A1     0
        0                γ2*k2/A2
        0                (1-γ2)k2/A3
        (1-γ1)k1/A4 0              ]

C = kc*[I(2) 0*I(2)] # Measure the first two tank levels
D = 0
G = ss(A,B,C,D)
nothing # hide
```
A PID controller with a filter is given by
```@example nonlinear
F = tf(1, [0.63, 1.12, 1])
Cpid = pid(0.26, 0.001, 15.9, form=:parallel)*F |> ss
nothing # hide
```
and to make the controller MIMO, we add a static pre-compensator that decouples the system at the the zero frequency.
```@example nonlinear
iG0 = dcgain(G)
iG0 ./= maximum(abs, iG0)
C = (Cpid .* I(2)) * iG0 
nothing # hide
```
The pumps (there are two of them) that service the tanks can only add liquid to the tanks, not remove liquid. The pump is thus saturated from below at 0, and from above at the maximum pump capacity 0.4. 
```@example nonlinear
using ControlSystemsBase: offset
umin = [0.0, 0.0]
umax = [0.4, 0.4]

yr    = G.C*xr  # Reference output
Gop   = offset(yr) * G * offset(-ur) # Make the plant operate in Δ-coordinates 
C_sat = saturation(umin, umax) * C   # while the controller and the saturation operate in the original coordinates
```
We now simulate the closed-loop system, the initial state of the plant is adjusted with the operating point `x0-xr` since the plant operates in Δ-coordinates 
```@example nonlinear
x0 = [2, 1, 8, 3] # Initial tank levels
plot(
    plot(lsim(feedback(Gop*C_sat), yr, 0:1:3000, x0=[x0-xr; zeros(C.nx)]), layout=1, sp=1, title="Outputs", ylabel=""),
    plot(lsim(feedback(C_sat, Gop), yr, 0:1:3000, x0=[zeros(C.nx); x0-xr]), layout=1, sp=1, title="Control signals", ylabel="")
)
hline!([yr[1]], label="Reference", l=:dash, sp=1, c=1)
```


### Duffing oscillator
In this example, we'll model and control the nonlinear system
```math
\ddot x = -kx - k_3 x^3 - c \dot{x} + 10u
```
To do this, we first draw the block diagram
```
10u    ┌───┐
──────►│+  │   ┌───┐   ┌───┐
 ┌────►│-  │ ẍ │ 1 │ ẋ │ 1 │ x
 │ ┌──►│-  ├──►│ - ├┬─►│ - ├─┬──►
 │ │ ┌►│-  │   │ s ││  │ s │ │
 │ │ │ └───┘   └───┘│  └───┘ │
 │ │ │              │        │
 │ │ │   ┌───┐      │        │
 │ │ └───┤ c │◄─────┘        │
 │ │     └───┘               │
 │ │                         │
 │ │     ┌───┐               │
 │ └─────┤ k │◄──────────────┤
 │       └───┘               │
 │                           │
 │       ┌───┐   ┌───┐       │
 └───────┤ k³│◄──┤ x³│◄──────┘
         └───┘   └───┘
```

We see that the input ``u`` passes through the inner velocity loop before reaching the output ``x``, we can form this inner closed-loop transfer function using `feedback(1/s, c)`, i.e., close the loop over an integrator by ``-c``. This inner loop is then connected in series with another integrator an feedback loop is closed with ``k_3 x^3 + kx = `` `pos_loop_feedback` in the feedback path. Notice how we multiply the final system with 10 from the right to get the input gain correct, for nonlinear systems, `10*sys` and `sys*10` are not always equivalent!

```@example DUFFING
using ControlSystems, Plots
using ControlSystemsBase: nonlinearity
k  = 10
k3 = 2
c  = 1

s = tf("s")

cube = nonlinearity(x->x^3)
vel_loop = feedback(1/s, c)
pos_loop_feedback = (k3*cube + k)
duffing = feedback(vel_loop/s, pos_loop_feedback)*10

plot(step(duffing, 20), title="Duffing oscillator open-loop step response")
```

We now show how we can make use of the circle criterion to prove stability of the closed loop. The function `circle_criterion` below plots the Nyquist curve of the loop-transfer function and figures out the circle to avoid by finding sector bounds for the static nonlinearity ``f(x) = x^3``. We then choose a controller and check that it stays outside of the circle. To find the sector bounds, we choose a domain to evaluate the nonlinearity over. The function ``f(x) = x^3`` goes to infinity faster than any linear function, and the upper sector bound is thus ∞, but if we restrict the nonlinearity to a smaller domain, we get a finite sector bound:
```@example DUFFING
function circle_criterion(L::ControlSystemsBase.HammersteinWienerSystem, domain::Tuple; N=10000)
    fun = x->L.f[](x)/x
    x = range(domain[1], stop=domain[2], length=N)
    0 ∈ x && (x = filter(!=(0), x)) # We cannot divide by zero
    k1, k2 = extrema(fun, x)

    f1 = plot(L.f[], domain[1], domain[2], title="Nonlinearity", lab="f(x)", xlab="x")
    plot!(x, [k1.*x k2.*x], lab=["k1 = $(round(k1, sigdigits=2))" "k2 = $(round(k2, sigdigits=2))"], l=(:dash), legend=:bottomright)

    p1 = -1/k2 # Close to origin
    p2 = -1/k1 # Far from origin

    c = (p1 + p2)/2
    r = (p2 - p1)/2

    Lnominal = sminreal(ss(L.A, L.B1, L.C1, L.D11, L.P.timeevol))
    f2 = nyquistplot(Lnominal)
    if p2 < -1000 # Due to bug in plots
        vspan!([-1000, p1], fillalpha=0.7, c=:red, primary=false)
    else
        th = 0:0.01:2pi
        Cs,Ss = cos.(th), sin.(th)
        plot!(r.*Cs .+ c, r.*Ss, fill=true, fillalpha=0.7, c=:red, primary=false)
    end

    plot(f1,f2)
end


C = pid(2, 0, 1, form=:parallel)*tf(1, [0.01,1])
f1 = circle_criterion(duffing*C, (-1, 1))
plot!(sp=2, ylims=(-10, 3), xlims=(-5, 11))
f2 = plot(step(feedback(duffing, C), 8), plotx=true, plot_title="Controlled oscillator disturbance step response", layout=4)
plot(f1,f2, size=(1300,800))
```
Since we evaluated the nonlinearity over a small domain, we should convince ourselves that we indeed never risk leaving this domain. 

In the example above, the circle turns into a half plane since the lower sector bound is 0. The example below chooses another nonlinearity
```math
f(x) = x + \sin(x)
```
to get an actual circle in the Nyquist plane.
```@example DUFFING
wiggly = nonlinearity(x->x+sin(x)) # This function is a bit wiggly
vel_loop = feedback(1/s, c)
pos_loop_feedback = (k3*wiggly + k)
duffing = feedback(vel_loop/s, pos_loop_feedback)*10

C = pid(2, 5, 1, form=:parallel)*tf(1,[0.1, 1]) 
f1 = circle_criterion(duffing*C, (-2pi, 2pi))
plot!(sp=2, ylims=(-5, 2), xlims=(-2.1, 0.1))
f2 = plot(step(feedback(duffing, C), 8), plotx=true, plot_title="Controlled wiggly oscillator disturbance step response", layout=5)
plot(f1,f2, size=(1300,800))
```


## Limitations
- Remember, this functionality is experimental and subject to breakage.
- Currently only `Continuous` systems supported.
- No nonlinear root-finding is performed during simulation. This limits the kinds of systems that can be simulated somewhat, in particular, no algebraic loops are allowed. 
- A lot of functions that expect linear systems will not work for nonlinear systems (naturally).

## Possible future work
- Discrete-time support.
- Basic support for nonlinear analysis such as stability proof through the circle criterion etc. In particular, predefined nonlinear functions may specify sector bounds for the gain, required by the circle-criterion calculations.
- Additional nonlinear components, such as 
    - Integrator anti-windup
    - Friction models

## See also
More advanced nonlinear modeling is facilitated by [ModelingToolkit.jl](https://github.com/SciML/ModelingToolkit.jl/) (MTK) and [ModelingToolkitStandardLibrary.jl](https://docs.sciml.ai/ModelingToolkitStandardLibrary/stable/). The tutorials 
- [Modeling for control using ModelingToolkit](https://help.juliahub.com/juliasimcontrol/dev/examples/mtk_control/)
- [Disturbance modeling in ModelingToolkit](https://help.juliahub.com/juliasimcontrol/dev/examples/mtk_disturbance_modeling/)
- [Modal analysis of a series of masses and springs using MTK](https://help.juliahub.com/juliasimcontrol/dev/linear_analysis/#Example:-Modal-analysis-of-a-series-of-masses-and-springs)

show how to use these packages to model and simulate control systems.

## Docstrings

```@docs
ControlSystemsBase.nonlinearity
ControlSystemsBase.offset
ControlSystemsBase.saturation
ControlSystemsBase.ratelimit
ControlSystemsBase.deadzone
```