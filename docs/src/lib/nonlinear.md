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
i.e., as a linear-fractional transform (LFT) between a linear system ``P`` and a diagonal matrix with scalar non-linear functions $f$. This representation is identical to that used for delay systems, and is exposed to the user in a similar way as well. The main entry point is the function [`nonlinearity`](@ref) which takes a nonlinear function $f$ like so, `nonlinearity(f)`. This creates a primitive system containing only the nonlinearity, but which behaves like a standard `LTISystem` during algebraic operations. We illustrate its usage through a number of examples:

## Examples
### Control-signal saturation
To create a controller that saturates the output at ``\pm 0.7``, we call
```@example nonlinear
using ControlSystems, Plots
using ControlSystems: nonlinearity # This functionality is not exported due to the beta status

C    = pid(; kp=1, ki=0.1)                  # A standard PI controller
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

Since the saturating nonlinearity is common, we provide the constructor [`saturation`](@ref) that automatically forms the equivalent to `nonlinearity(x->clamp(x, -0.7, 0.7))` while at the same time making sure the function has a recognizable name when the system is printed
```@example nonlinear
using ControlSystems: saturation
saturation(0.7)
```

### Non-zero operating point
It's common to linearize nonlinear systems around some operating point. We may make use of the helper constructor [`offset`](@ref) to create affine functions at the inputs and outputs of the linearized system to, e.g.,
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
Cpid = ss(pid(;kp=0.26, ki=0.001, kd=15.9)*F)
nothing # hide
```
and to make the controller MIMO, we add a static pre-compensator that decouples the system at the the zero frequency.
```@example nonlinear
iG0 = dcgain(G)
iG0 ./= maximum(abs, iG0)
C = Cpid*iG0
nothing # hide
```
The pumps (there are two of them) that service the tanks can only add liquid to the tanks, not remove liquid. The pump is thus saturated from below at 0, and from above at the maximum pump capacity 0.4. 
```@example nonlinear
using ControlSystems: offset
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

## Limitations
- Remember, this functionality is experimental and subject to breakage.
- Currently only `Continuous` systems supported. Discrete-time systems will come soon.
- No nonlinear root-finding is performed during simulation. This limits the kinds of systems that can be simulated somewhat, in particular, no algebraic loops are allowed. 
- A lot of functions that expect linear systems will not work for nonlinear systems (naturally).

## Future work
- Discrete-time support.
- Basic support for nonlinear analysis such as stability proof through the circle criterion etc. In particular, predefined nonlinear functions may specify sector bounds for the gain, required by the circle-criterion calculations.
- Additional nonlinear components, such as 
    - Dead zone
    - Integrator anti-windup
    - Friction models?

## Docstrings

```@docs
ControlSystems.nonlinearity
ControlSystems.offset
ControlSystems.saturation
```