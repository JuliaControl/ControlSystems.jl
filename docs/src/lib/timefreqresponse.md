```@index
Pages = ["timefreqresponse.md"]
```

# Time and Frequency response

Any `TransferFunction` can be evaluated at a point using
`F(s)`, `F(omega, true)`, `F(z, false)`

- `F(s)` evaluates the continuous-time transfer function `F` at `s`.
- `F(omega,true)` evaluates the discrete-time transfer function `F` at `exp(i*Ts*omega)`
- `F(z,false)` evaluates the discrete-time transfer function `F` at `z`

```@autodocs
Modules = [ControlSystems]
Pages   = ["timeresp.jl", "freqresp.jl", "simulators.jl"]
Order   = [:function, :type]
Private = false
```
