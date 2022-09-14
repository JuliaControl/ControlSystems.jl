```@index
Pages = ["plotting.md"]
```

!!! note "Using Plots"
    All plotting requires the user to manually load the Plots.jl library, e.g., by calling `using Plots`.

!!! note "Time-domain responses"
    There are no special functions to plot time-domain results, such as step and impulse responses, instead, simply call `plot` on the result structure ([`ControlSystemsBase.SimResult`](@ref)) returned by [`lsim`](@ref), [`step`](@ref), [`impulse`](@ref) etc.

# Plotting functions

```@autodocs
Modules = [ControlSystems, ControlSystemsBase]
Pages   = ["lib/ControlSystemsBase/src/plotting.jl"]
Order   = [:function]
Private = false
```

## Examples

### Bode plot

![bode](https://github.com/JuliaControl/ControlExamplePlots.jl/blob/master/src/figures/bode.png?raw=true)

### Sigma plot

![sigma](https://github.com/JuliaControl/ControlExamplePlots.jl/blob/master/src/figures/sigma.png?raw=true)

### Margin

![margin](https://github.com/JuliaControl/ControlExamplePlots.jl/blob/master/src/figures/margin.png?raw=true)

### Gangoffour plot

![gangoffour](https://github.com/JuliaControl/ControlExamplePlots.jl/blob/master/src/figures/gangoffour.png?raw=true)

### Nyquist plot

![nyquist](https://github.com/JuliaControl/ControlExamplePlots.jl/blob/master/src/figures/nyquist.png?raw=true)

### Nichols plot

![nichols](https://github.com/JuliaControl/ControlExamplePlots.jl/blob/master/src/figures/nichols.png?raw=true)

### Pole-zero plot

![pzmap](https://github.com/JuliaControl/ControlExamplePlots.jl/blob/master/src/figures/pzmap.png?raw=true)

### Rlocus plot

![rlocus](https://github.com/JuliaControl/ControlExamplePlots.jl/blob/master/src/figures/rlocus.png?raw=true)

### Lsim response plot

![lsim](https://github.com/JuliaControl/ControlExamplePlots.jl/blob/master/src/figures/lsim.png?raw=true)

### Impulse response plot

![impulse](https://github.com/JuliaControl/ControlExamplePlots.jl/blob/master/src/figures/impulse.png?raw=true)

### Step response plot

![step](https://github.com/JuliaControl/ControlExamplePlots.jl/blob/master/src/figures/step.png?raw=true)

