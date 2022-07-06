# Creating Systems
## Creating Transfer Functions
```@meta
DocTestSetup = quote
    using ControlSystems
end
```

### tf - Rational Representation
The syntax for creating a transfer function is
```julia
tf(num, den)     # Continuous-time system
tf(num, den, Ts) # Discrete-time system
```
where `num` and `den` are the polynomial coefficients of the numerator and denominator of the polynomial and `Ts`, if provided, is the sample time for a discrete-time system.
#### Example:
```jldoctest
tf([1.0],[1,2,1])

# output

TransferFunction{Continuous, ControlSystems.SisoRational{Float64}}
        1.0
-------------------
1.0s^2 + 2.0s + 1.0

Continuous-time transfer function model
```

The transfer functions created using this method will be of type `TransferFunction{SisoRational}`.

### zpk - Pole-Zero-Gain Representation
Sometimes it's better to represent the transfer function by its poles, zeros and gain, this can be done using
```julia
zpk(zeros, poles, gain)     # Continuous-time system
zpk(zeros, poles, gain, Ts) # Discrete-time system
```
where `zeros` and `poles` are `Vectors` of the zeros and poles for the system and `gain` is a gain coefficient.
#### Example
```jldoctest
zpk([-1.0,1], [-5, -10], 2)

# output

TransferFunction{Continuous, ControlSystems.SisoZpk{Float64, Float64}}
   (1.0s + 1.0)(1.0s - 1.0)
2.0-------------------------
   (1.0s + 5.0)(1.0s + 10.0)

Continuous-time transfer function model
```

The transfer functions created using this method will be of type `TransferFunction{SisoZpk}`.


## ss - Creating State-Space Systems
A state-space system is created using
```julia
ss(A,B,C,D)    # Continuous-time system
ss(A,B,C,D,Ts) # Discrete-time system
```
and they behave similarily to transfer functions. State-space systems with heterogeneous matrix types are also available, which can be used to create systems with static or sized matrices, e.g.,
```jldoctest HSS; output=false
using StaticArrays
import ControlSystems.HeteroStateSpace
@inline to_static(a::Number) = a
@inline to_static(a::AbstractArray) = SMatrix{size(a)...}(a)
@inline to_sized(a::Number) = a
@inline to_sized(a::AbstractArray) = SizedArray{Tuple{size(a)...}}(a)
function HeteroStateSpace(A,B,C,D,Ts=0,f::F=to_static) where F
    HeteroStateSpace(f(A),f(B),f(C),f(D),Ts)
end
@inline HeteroStateSpace(s,f) = HeteroStateSpace(s.A,s.B,s.C,s.D,s.timeevol,f)
ControlSystems._string_mat_with_headers(a::SizedArray) = ControlSystems._string_mat_with_headers(Matrix(a)); # Overload for printing purposes

# output

```
Notice the different matrix types used
```jldoctest HSS
julia> sys = ss([-5 0; 0 -5],[2; 2],[3 3],[0])
StateSpace{Continuous, Int64}
A =
 -5   0
  0  -5
B =
 2
 2
C =
 3  3
D =
 0

Continuous-time state-space model

julia> HeteroStateSpace(sys, to_static)
HeteroStateSpace{Continuous, SMatrix{2, 2, Int64, 4}, SMatrix{2, 1, Int64, 2}, SMatrix{1, 2, Int64, 2}, SMatrix{1, 1, Int64, 1}}
A =
 -5   0
  0  -5
B =
 2
 2
C =
 3  3
D =
 0

Continuous-time state-space model

julia> HeteroStateSpace(sys, to_sized)
HeteroStateSpace{Continuous, SizedMatrix{2, 2, Int64, 2, Matrix{Int64}}, SizedMatrix{2, 1, Int64, 2, Matrix{Int64}}, SizedMatrix{1, 2, Int64, 2, Matrix{Int64}}, SizedMatrix{1, 1, Int64, 2, Matrix{Int64}}}
A =
 -5   0
  0  -5
B =
 2
 2
C =
 3  3
D =
 0

Continuous-time state-space model
```

## Converting between types
It is sometime useful to convert one representation to another, this is possible using the same functions, for example
```jldoctest
tf(zpk([-1], [1], 2, 0.1))

# output

TransferFunction{Discrete{Float64}, ControlSystems.SisoRational{Int64}}
2z + 2
------
z - 1

Sample Time: 0.1 (seconds)
Discrete-time transfer function model
```

## Creating Delay Systems
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

## Creating Nonlinear Systems
See [Nonlinear functionality](@ref).

## Simplifying systems
A statespace system with a non-minimal realization, or a transfer function with overlapping zeros and poles, may be simplified using the function [`minreal`](@ref). Systems that are structurally singular, i.e., that contains outputs that can not be reached from the inputs based on analysis of the structure of the zeros in the system matrices only, can be simplified with the function [`sminreal`](@ref).

Examples:
```@repl
using ControlSystems
G = tf([1, 1], [1, 1])
minreal(G) # Performs pole-zero cancellation

P = tf(1, [1, 1]) |> ss
G = P / (1 + P) # this creates a non-minimal realization, use feedback(P) instead
feedback(P) # Creates a minimal realization directly
Gmin = minreal(G) # this simplifies the realization to a minimal realization
norm(Gmin - feedback(P), Inf) # No difference
bodeplot([G, Gmin, feedback(P)]) # They are all identical
```