```@index
Pages = ["timefreqresponse.md"]
```

# Time and Frequency response

Frequency responses are calculated using [`freqresp`](@ref), [`bode`](@ref) and [`nyquist`](@ref).

Any [`TransferFunction`](@ref) can be evaluated at a point using
`F(s)`, `F(omega, true)`, `F(z, false)`

- `F(s)` evaluates the continuous-time transfer function `F` at `s`.
- `F(omega,true)` evaluates the discrete-time transfer function `F` at `exp(i*Ts*omega)`
- `F(z,false)` evaluates the discrete-time transfer function `F` at `z`

Simulation with arbitrary inputs is primarily handled by the function [`lsim`](@ref), with [`step`](@ref) and [`impulse`](@ref) serving as convenience functions to simulate responses to particular inputs.

For more extensive nonlinear simulation capabilities, see the notes on ModelingToolkit and DifferentialEquations under [The wider Julia ecosystem for control](@ref).

```@autodocs
Modules = [ControlSystems, ControlSystemsBase]
Pages   = [libpath*"/timeresp.jl", libpath*"/result_types.jl", libpath*"/freqresp.jl", "simulators.jl"]
Order   = [:function, :type]
Private = false
```
