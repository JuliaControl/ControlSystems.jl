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
where `num` and `den` are the polinomial coefficients of the numerator and denominator of the polynomial and `Ts` is the sample time.
### Example:
```julia
tf([1.0],[1,2,1])

# output

TransferFunction{ControlSystems.SisoRational{Float64}}
         1.0
---------------------
1.0*s^2 + 2.0*s + 1.0


Continuous-time transfer function model
```

The transfer functions created using this method will be of type `TransferFunction{SisoRational}`.

## zpk - Pole-Zero-Gain Representation
Sometimes it's better to represent the transferfunction by its poles, zeros and gain, this can be done using
```julia
zpk(zeros, poles, gain, Ts=0)
```
where `zeros` and `poles` are `Vectors` of the zeros and poles for the system and `gain` is a gain coefficient.
### Example
```julia
zpk([-1.0,1], [-5, -10], 2)

# output

TransferFunction{ControlSystems.SisoZpk{Float64,Float64}}
   (1.0*s + 1.0)(1.0*s - 1.0)
2.0---------------------------
   (1.0*s + 5.0)(1.0*s + 10.0)

Continuous-time transfer function model
```

The transfer functions created using this method will be of type `TransferFunction{SisoZpk}`.

## Converting between types
It is sometime useful to convert one representation to another, this is possible using the same functions, for example
```julia
tf(zpk([-1], [1], 2, 0.1))

# output

TransferFunction{ControlSystems.SisoRational{Int64}}
2*z + 2
-------
1z - 1

Sample Time: 0.1 (seconds)
Discrete-time transfer function model
```
