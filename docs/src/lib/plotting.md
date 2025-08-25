```@index
Pages = ["plotting.md"]
```

!!! note "Using Plots"
    All plotting requires the user to manually load the Plots.jl library, e.g., by calling `using Plots`.

!!! note "Time-domain responses"
    There are no special functions to plot time-domain results, such as step and impulse responses, instead, simply call `plot` on the result structure (`ControlSystemsBase.SimResult`) returned by [`lsim`](@ref), [`step`](@ref), [`impulse`](@ref) etc.

## Makie support
!!! danger "Experimental"
    
    The support for plotting with Makie is currently experimental and at any time subject to breaking changes or removal **not** respecting semantic versioning.

ControlSystemsBase provides experimental support for plotting with [Makie.jl](https://docs.makie.org/) through the `CSMakie` module. This support is loaded automatically when you load a Makie backend (GLMakie, CairoMakie, or WGLMakie).

### Usage

```julia
using ControlSystemsBase, GLMakie  # or CairoMakie, WGLMakie

# Create a system
P = tf([1], [1, 2, 1])

# Use CSMakie plotting functions
CSMakie.bodeplot(P)
CSMakie.nyquistplot(P)
CSMakie.pzmap(P)
# ... and more

# Direct plotting of simulation results
res = step(P, 10)
plot(res)  # Creates a figure with time-domain response

si = stepinfo(res)
plot(si)  # Visualizes step response characteristics
```

### Available functions

The `CSMakie` module provides Makie implementations of the following plotting functions:

- `CSMakie.bodeplot` - Bode magnitude and phase plots
- `CSMakie.nyquistplot` - Nyquist plots with optional M and Mt circles
- `CSMakie.sigmaplot` - Singular value plots
- `CSMakie.marginplot` - Gain and phase margin plots
- `CSMakie.pzmap` - Pole-zero maps
- `CSMakie.nicholsplot` - Nichols charts
- `CSMakie.rgaplot` - Relative gain array plots
- `CSMakie.rlocusplot` - Root locus plots
- `CSMakie.leadlinkcurve` - Lead-link design curves

Additionally, `SimResult` and `StepInfo` types can be plotted directly using Makie's `plot` function.

# Plotting functions

```@autodocs
Modules = [ControlSystems, ControlSystemsBase]
Pages   = [libpath*"/plotting.jl"]
Order   = [:function]
Private = false
```
```@docs
ControlSystemsBase.rlocusplot
```
- To plot simulation results such as step and impulse responses, use `plot(::SimResult)`, see also [`lsim`](@ref).

## Examples

### Bode plot

![bode](https://github.com/JuliaControl/ControlExamplePlots.jl/blob/master/src/figures/bode.png?raw=true)
```julia
tf1 = tf([1],[1,1])
tf2 = tf([1/5,2],[1,1,1])
sys = [tf1 tf2]
ws = exp10.(range(-2,stop=2,length=200))
bodeplot(sys, ws)
```
### Sigma plot

![sigma](https://github.com/JuliaControl/ControlExamplePlots.jl/blob/master/src/figures/sigma.png?raw=true)
```julia
sys = ss([-1 2; 0 1], [1 0; 1 1], [1 0; 0 1], [0.1 0; 0 -0.2])
sigmaplot(sys)
```
### Margin

![margin](https://github.com/JuliaControl/ControlExamplePlots.jl/blob/master/src/figures/margin.png?raw=true)
```julia
tf1 = tf([1],[1,1])
tf2 = tf([1/5,2],[1,1,1])
ws = exp10.(range(-2,stop=2,length=200))
marginplot([tf1, tf2], ws)
```
### Gangoffour plot

![gangoffour](https://github.com/JuliaControl/ControlExamplePlots.jl/blob/master/src/figures/gangoffour.png?raw=true)
```julia
tf1 = tf([1.0],[1,1])
gangoffourplot(tf1, [tf(1), tf(5)])
```
### Nyquist plot

![nyquist](https://github.com/JuliaControl/ControlExamplePlots.jl/blob/master/src/figures/nyquist.png?raw=true)
```julia
sys = ss([-1 2; 0 1], [1 0; 1 1], [1 0; 0 1], [0.1 0; 0 -0.2])
ws = exp10.(range(-2,stop=2,length=200))
nyquistplot(sys, ws, Ms_circles=1.2, Mt_circles=1.2)
```
### Nichols plot

![nichols](https://github.com/JuliaControl/ControlExamplePlots.jl/blob/master/src/figures/nichols.png?raw=true)
```julia
tf1 = tf([1],[1,1])
ws = exp10.(range(-2,stop=2,length=200))
nicholsplot(tf1,ws)
```
### Pole-zero plot

![pzmap](https://github.com/JuliaControl/ControlExamplePlots.jl/blob/master/src/figures/pzmap.png?raw=true)
```julia
tf2 = tf([1/5,2],[1,1,1])
pzmap(c2d(tf2, 0.1))
```
### Rlocus plot

![rlocus](https://github.com/JuliaControl/ControlExamplePlots.jl/blob/master/src/figures/rlocus.png?raw=true)

### Lsim response plot

![lsim](https://github.com/JuliaControl/ControlExamplePlots.jl/blob/master/src/figures/lsim.png?raw=true)
```julia
sys = ss([-1 2; 0 1], [1 0; 1 1], [1 0; 0 1], [0.1 0; 0 -0.2])
sysd = c2d(sys, 0.01)
L = lqr(sysd, [1 0; 0 1], [1 0; 0 1])
ts = 0:0.01:5
plot(lsim(sysd, (x,i)->-L*x, ts; x0=[1;2]), plotu=true)
```
### Impulse response plot

![impulse](https://github.com/JuliaControl/ControlExamplePlots.jl/blob/master/src/figures/impulse.png?raw=true)
```julia
tf1 = tf([1],[1,1])
tf2 = tf([1/5,2],[1,1,1])
sys = [tf1 tf2]
sysd = c2d(ss(sys), 0.01)
plot(impulse(sysd, 5), l=:blue)
```
### Step response plot

![step](https://github.com/JuliaControl/ControlExamplePlots.jl/blob/master/src/figures/step.png?raw=true)
```julia
tf1 = tf([1],[1,1])
tf2 = tf([1/5,2],[1,1,1])
sys = [tf1 tf2]
sysd = c2d(ss(sys), 0.01)
res = step(sysd, 5)
plot(res, l=(:dash, 4))
# plot!(stepinfo(step(sysd[1,1], 5))) # adds extra info to the plot
```
