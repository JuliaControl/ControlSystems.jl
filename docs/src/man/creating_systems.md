# Creating Transfer Functions
```@meta
DocTestSetup = quote
    using ControlSystems
end
```

## tf - Rational Representation
The syntax for creating a transfer function is
```julia
tf(num, den, Ts=0)
```
where `num` and `den` are the polynomial coefficients of the numerator and denominator of the polynomial and `Ts` is the sample time.
### Example:
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

## zpk - Pole-Zero-Gain Representation
Sometimes it's better to represent the transfer function by its poles, zeros and gain, this can be done using
```julia
zpk(zeros, poles, gain, Ts=0)
```
where `zeros` and `poles` are `Vectors` of the zeros and poles for the system and `gain` is a gain coefficient.
### Example
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


# Creating State-Space Systems
A state-space system is created using
```julia
ss(A,B,C,D,Ts=0)
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
HeteroStateSpace{Continuous, SizedMatrix{2, 2, Int64, 2}, SizedMatrix{2, 1, Int64, 2}, SizedMatrix{1, 2, Int64, 2}, SizedMatrix{1, 1, Int64, 2}}
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
