# Creating Transfer Functions
    {meta}
    DocTestSetup = quote
        using ControlSystems
    end

## tf - Rational Representation
The syntax for creating a transfer function is
```julia
tf(num, den, Ts=0)
```
where `num` and `den` are the polinomial coefficients of the numerator and denominator of the polynomial and `Ts` is the sample time.
### Example:
```julia
tf([1],[1,2,1])

# output

TransferFunction:
      1.0
----------------
s^2 + 2.0s + 1.0

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
zpk([-1,1], [-5, -10], 2)

# output

TransferFunction:
   (s - 1.0)(s + 1.0)
2.0-------------------
   (s + 10.0)(s + 5.0)

Continuous-time transfer function model
```

The transfer functions created using this method will be of type `TransferFunction{SisoZpk}`.

## tfa - Generalized Representation
If you want to work with transfer functions that are not rational functions, it is possible to use the `tfa` representation
```julia
tfa(str::String), tfa(str::Expr)
```
This function will either convert `str` to an expression or directly accept an `Expr` and create a transfer function.
### Example:
```julia
tfa("1/((s+1)*exp(-sqrt(s)))")

## output

TransferFunction:
1/((s+1)*exp(-sqrt(s)))

Continuous-time transfer function model
```
The transfer functions created using this method will be of type `TransferFunction{SisoGeneralized}`.
This type will work with some functions like `bodeplot, stepplot` but not others ,like `poles`.

## Converting between types
It is sometime useful to convert one representation to another, this is possible using the same functions, for example
```julia
tf(zpk([-1], [1], 2, 0.1))

# output

TransferFunction:
2.0z + 2.0
----------
 z - 1.0

Sample Time: 0.1 (seconds)
Discrete-time transfer function model
```
