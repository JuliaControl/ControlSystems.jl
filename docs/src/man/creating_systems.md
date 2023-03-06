# Creating Systems
This page illustrates how to create system models such as transfer functions and statespace models. This topic is also treated in the introductory video below:

```@raw html
<iframe style="height: 315px; width: 560px" src="https://www.youtube.com/embed/ksrEyMNX_BY" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>
```

## Transfer Functions
```@meta
DocTestSetup = quote
    using ControlSystems
end
```

### tf - Rational Representation
The syntax for creating a transfer function is [`tf`](@ref)
```julia
tf(num, den)     # Continuous-time system
tf(num, den, Ts) # Discrete-time system
```
where `num` and `den` are the polynomial coefficients of the numerator and denominator of the polynomial and `Ts`, if provided, is the sample time for a discrete-time system.
#### Example:
```jldoctest
tf([1.0],[1,2,1])

# output

TransferFunction{Continuous, ControlSystemsBase.SisoRational{Float64}}
        1.0
-------------------
1.0s^2 + 2.0s + 1.0

Continuous-time transfer function model
```

The transfer functions created using this method will be of type `TransferFunction{SisoRational}`.

### zpk - Pole-Zero-Gain Representation
Sometimes it's better to represent the transfer function by its poles, zeros and gain, this can be done using the function [`zpk`](@ref)
```julia
zpk(zeros, poles, gain)     # Continuous-time system
zpk(zeros, poles, gain, Ts) # Discrete-time system
```
where `zeros` and `poles` are `Vectors` of the zeros and poles for the system and `gain` is a gain coefficient.
#### Example
```jldoctest
zpk([-1.0,1], [-5, -10], 2)

# output

TransferFunction{Continuous, ControlSystemsBase.SisoZpk{Float64, Float64}}
   (1.0s + 1.0)(1.0s - 1.0)
2.0-------------------------
   (1.0s + 5.0)(1.0s + 10.0)

Continuous-time transfer function model
```

The transfer functions created using this method will be of type `TransferFunction{SisoZpk}`.


## State-Space Systems
A state-space system is created using
```julia
ss(A,B,C,D)    # Continuous-time system
ss(A,B,C,D,Ts) # Discrete-time system
```
and they behave similarly to transfer functions.

The [`ss`](@ref) constructor allows you to
- Pass `0` instead of a ``D`` matrix, and an appropriately sized zero matrix is created automatically.
- Pass `I` instead of a ``C`` matrix, and an appropriately sized identity matrix is created automatically. The `UniformScaling` operator `I` lives in the `LinearAlgebra` standard library which must be loaded first.

State-space systems with heterogeneous matrix types are also available, which can be used to create systems with static or sized matrices, e.g.,
```@example HSS
using ControlSystemsBase, StaticArrays
sys = ss([-5 0; 0 -5],[2; 2],[3 3],[0])
HeteroStateSpace(sys, to_sized)
HeteroStateSpace(sys, to_static)
```
Notice the different matrix types used.

To associate **names** with states, inputs and outputs, see [`named_ss`](https://juliacontrol.github.io/RobustAndOptimalControl.jl/dev/#Named-systems) from RobustAndOptimalControl.jl.


## Converting between types
It is sometime useful to convert one representation to another, this is possible using the constructors `tf, zpk, ss`, for example
```jldoctest
tf(zpk([-1], [1], 2, 0.1))

# output

TransferFunction{Discrete{Float64}, ControlSystemsBase.SisoRational{Int64}}
2z + 2
------
z - 1

Sample Time: 0.1 (seconds)
Discrete-time transfer function model
```

## Delay Systems
The constructor [`delay`](@ref) creates a pure delay, which may be connected to a system by multiplication:
```julia
delay(1.2)               # Pure delay or 1.2s
tf(1, [1, 1])*delay(1.2) # Input delay
delay(1.2)*tf(1, [1, 1]) # Output delay
```

Delayed systems can also be created using
```julia
s = tf("s")
L = 1.2 # Delay time
tf(1, [1, 1]) * exp(-L*s)
```

Padé approximations of delays can be created using [`pade`](@ref).

A tutorial on delay systems is available here:
```@raw html
<iframe style="height: 315px; width: 560px" src="https://www.youtube.com/embed/ER8_oHU2vZs" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>
```


## Nonlinear Systems
See [Nonlinear functionality](@ref).

## Simplifying systems
A statespace system with a non-minimal realization, or a transfer function with overlapping zeros and poles, may be simplified using the function [`minreal`](@ref). Systems that are structurally singular, i.e., that contains outputs that can not be reached from the inputs based on analysis of the structure of the zeros in the system matrices only, can be simplified with the function [`sminreal`](@ref).

Examples:
```@repl
using ControlSystemsBase
G = tf([1, 1], [1, 1])
minreal(G) # Performs pole-zero cancellation

P = tf(1, [1, 1]) |> ss
G = P / (1 + P) # this creates a non-minimal realization, use feedback(P) instead
feedback(P) # Creates a minimal realization directly
Gmin = minreal(G) # this simplifies the realization to a minimal realization
norm(Gmin - feedback(P), Inf) # No difference
bodeplot([G, Gmin, feedback(P)]) # They are all identical
```

## Multiplying systems
Two systems can be connected in series by multiplication
```@example MIMO
using ControlSystemsBase
P1 = ss(-1,1,1,0)
P2 = ss(-2,1,1,0)
P2*P1
```
If the input dimension of `P2` does not match the output dimension of `P1`, an error is thrown. If one of the systems is SISO and the other is MIMO, broadcasted multiplication will expand the SISO system to match the input or output dimension of the MIMO system, e.g.,
```@example MIMO
Pmimo = ssrand(2,2,1)
Psiso = ss(-2,1,1,0)
# Psiso * Pmimo # error
Psiso .* Pmimo ≈ [Psiso 0; 0 Psiso] * Pmimo # Broadcasted multiplication expands SISO into diagonal system
```

Broadcasted multiplication between a system and an array is only allowed for diagonal arrays
```@example MIMO
using LinearAlgebra
Psiso .* I(2)
```


## MIMO systems and arrays of systems
Concatenation of systems creates MIMO systems, which is different from an array of systems. For example
```@example MIMO
using ControlSystemsBase
P = ss(-1,1,1,0)
P_MIMO = [P 2P]
```
is a 1×2 MISO system, not a 1×2 array.

### From SISO to MIMO
SISO systems do not multiply MIMO systems directly, i.e.,
```@example MIMO
using Test
siso = ss(-1,1,1,0)
mimo = ssrand(2,2,2)
@test_throws DimensionMismatch siso * mimo
```

To multiply `siso` with each output channel of `mimo` in the example above, use broadcasting:
```@example MIMO
siso .* mimo
```
This is equivalent to first expanding the SISO system into a diagonal system
```@example MIMO
using LinearAlgebra
(siso .* I(2)) * mimo
```



### Converting an array of systems to a MIMO system
Diagonal MIMO systems can be created from a vector of systems using [`append`](@ref)
```@example MIMO
P1 = ssrand(1,1,1)
P2 = ssrand(1,1,1)
append(P1, P2)
```
More general arrays of systems can be converted to a MIMO system using [`array2mimo`](@ref).
```@example MIMO
sys_array = fill(P, 2, 2) # Creates an array of systems
mimo_sys = array2mimo(sys_array)
```

### Converting MIMO system to an array of systems
This conversion is not explicitly supported, but is easy enough to accomplish with standard Julia code, for example:
```@example MIMO
P = ssrand(2,3,1) # A random 2×3 MIMO system
sys_array = getindex.(Ref(P), 1:P.ny, (1:P.nu)')
```

### Creating arrays with different types of systems
When calling `hcat/vcat`, Julia automatically tries to promote the types to the smallest common supertype, this means that creating an array with one continuous and one discrete-time system fails
```@example MIMO
P_cont = ssrand(2,3,1) 
P_disc = ssrand(2,3,1, Ts=1)
@test_throws ErrorException [P_cont, P_disc] # ERROR: Sampling time mismatch
```
You can explicitly tell Julia that you want a particular supertype, e.g,
```@example MIMO
StateSpace[P_cont, P_disc]
```
The type `StateSpace` is abstract, since the type parameters are not specified.

## Demo systems
The module `ControlSystemsBase.DemoSystems` contains a number of demo systems demonstrating different kinds of dynamics.


## From block diagrams to code
This section lists a number of block diagrams, and indicates the corresponding transfer functions and how they are built in code.

The function `feedback(G1, G2)` can be thought of like this: the first argument `G1` is the system that appears directly between the input and the output (the *forward path*), while the second argument `G2` (defaults to 1 if omitted) contains all other systems that appear in the closed loop (the *feedback path*). The feedback is assumed to be negative, unless the argument `pos_feedback = true` is passed (`lft` is an exception, which due to convention defaults to positive feedback). This means that `feedback(G, 1)` results in unit negative feedback, while `feedback(G, -1)` or `feedback(G, 1, pos_feedback = true)` results in unit positive feedback.

---
Closed-loop system from reference to output
```
r   ┌─────┐     ┌─────┐
───►│     │  u  │     │ y
    │  C  ├────►│  P  ├─┬─►
 -┌►│     │     │     │ │
  │ └─────┘     └─────┘ │
  │                     │
  └─────────────────────┘
```
```math
Y = \dfrac{PC}{I+PC}R
```

Code: `feedback(P*C)` or equivalently `comp_sensitivity(P, C)`. Here, the system ``PC`` appears directly between the input ``r`` and the output ``y``, and the feedback loop is negative identity. We thus call `feedback(P*C) = feedback(P*C, 1)`

---

```
d     ┌───┐   y
───+─►│ P ├─┬───►
  -▲  └───┘ │
   │        │
   │  ┌───┐ │
   └──┤ C │◄┘
      └───┘
```
```math
Y = \dfrac{P}{I+PC}D = PSD
```

Code: `feedback(P, C)` or equivalently `G_PS(P, C)`. Here, only ``P`` appears directly between ``d`` and ``y``, while ``C`` appears first in the feedback loop. We thus call `feedback(P, C)`

---
Sensitivity function at plant input

```
d    e┌───┐   
───+─►│ P ├─┬───►
  -▲  └───┘ │
   │        │
   │  ┌───┐ │
   └──┤ C │◄┘
      └───┘
```
```math
E = \dfrac{1}{I+CP}D = SD
```

Code: `feedback(1, C*P)` or equivalently `input_sensitivity(P, C)`. Here, there are no systems directly between the input and the output, we thus call `feedback(1, C*P)`. Note the order in `C*P`, which is important for MIMO systems. This computes the sensitivity function at the *plant input*. It's more common to analyze the sensitivity function at the plant output, illustrated below (for SISO systems they are equivalent).

---

Sensitivity function at plant output

```
      ┌───┐   
───+─►│ P ├─+◄── e
  -▲  └───┘ │
   │        │y
   │  ┌───┐ │
   └──┤ C │◄┘
      └───┘
```
```math
Y = \dfrac{1}{I+PC}E = SE
```

Code: `feedback(1, P*C)` or equivalently `output_sensitivity(P, C)`. Note the reverse order in ``PC`` compared to the input sensitivity function above.


---

Reference ``r`` and input disturbance ``d`` to output ``y`` and control signal ``u``. This example forms the transfer function matrix with ``r`` and ``d`` as inputs, and ``y`` and ``u`` as outputs.
```
              d
     ┌─────┐  │  ┌─────┐
r    │     │u ▼  │     │ y
──+─►│  C  ├──+─►│  P  ├─┬─►
  ▲  │     │     │     │ │
 -│  └─────┘     └─────┘ │
  │                      │
  └──────────────────────┘
```

```math
\begin{bmatrix}
y \\ u
\end{bmatrix} = 
\begin{bmatrix}
\dfrac{PC}{I + PC} & \dfrac{C}{I + PC} \\
\dfrac{P}{I + PC} & \dfrac{-PC}{I + PC}
\end{bmatrix}
\begin{bmatrix}
r \\ d
\end{bmatrix}
```

Code: `feedback(C, P, W2=:, Z2=:, Zperm=[(1:P.ny).+P.nu; 1:P.nu]) # y,u from r,d`.
Here, we have reversed the order of `P` and `C` to get the correct sign of the control signal. We also make use of the keyword arguments `W2` and `Z2` to specify that we want to include the inputs and outputs of `P` as external inputs and outputs, and `Zperm` to specify the order of the outputs (``y`` before ``u``).

---

Linear fractional transformation

```
     ┌─────────┐
z◄───┤         │◄────w
     │    P    │
y┌───┤         │◄───┐u
 │   └─────────┘    │
 │                  │
 │      ┌───┐       │
 │      │   │       │
 └─────►│ K ├───────┘
        │   │
        └───┘
```
```math
Z = \operatorname{lft}{(P, K)} W
```

Code: `lft(P, K)`

---

```
      z1          z2
      ▲  ┌─────┐  ▲      ┌─────┐
      │  │     │  │      │     │
w1──+─┴─►│  C  ├──┴───+─►│  P  ├─┐
    │    │     │      │  │     │ │
    │    └─────┘      │  └─────┘ │
    │                 w2         │
    └────────────────────────────┘
```

The transfer function from ``w_1, w_2`` to ``z_1, z_2`` contains all the transfer functions that are commonly called "gang of four" (see also [`gangoffour`](@ref)).

```math
\begin{bmatrix}
z_1 \\ z_2
\end{bmatrix} = 
\begin{bmatrix}
I \\ C
\end{bmatrix} (I + PC)^{-1} \begin{bmatrix}
I & P
\end{bmatrix}
\begin{bmatrix}
w_1 \\ w_2
\end{bmatrix}
```
Code: This function requires the package [RobustAndOptimalControl.jl](https://juliacontrol.github.io/RobustAndOptimalControl.jl/dev/).
```julia
RobustAndOptimalControl.extended_gangoffour(P, C, pos=true)
# For SISO P
S  = G[1, 1]
PS = G[1, 2]
CS = G[2, 1]
T  = G[2, 2]

# For MIMO P
S  = G[1:P.ny,     1:P.nu]
PS = G[1:P.ny,     P.nu+1:end]
CS = G[P.ny+1:end, 1:P.nu]
T  = G[P.ny+1:end, P.nu+1:end]
```

See also
- [`output_sensitivity`](@ref)
- [`input_sensitivity`](@ref)
- [`output_comp_sensitivity`](@ref)
- [`input_comp_sensitivity`](@ref)
- [`G_PS`](@ref)
- [`G_CS`](@ref)
- [`gangoffour`](@ref))
- [`gangoffourplot`](@ref))

---

This diagram is more complicated and forms several connections, including both feedforward and feedback connections. A code file that goes through how to form such complicated connections using named signals is linked below. This example uses the package [RobustAndOptimalControl.jl](https://juliacontrol.github.io/RobustAndOptimalControl.jl/dev/).

```
                 yF
              ┌────────────────────────────────┐
              │                                │
    ┌───────┐ │  ┌───────┐ yR   ┌─────────┐    │    ┌───────┐
uF  │       │ │  │       ├──────►         │ yC │  uP│       │    yP
────►   F   ├─┴──►   R   │      │    C    ├────+────►   P   ├────┬────►
    │       │    │       │   ┌──►         │         │       │    │
    └───────┘    └───────┘   │- └─────────┘         └───────┘    │
                             │                                   │
                             └───────────────────────────────────┘
```

See code example [complicated_feedback.jl](https://github.com/JuliaControl/RobustAndOptimalControl.jl/blob/master/examples/complicated_feedback.jl).