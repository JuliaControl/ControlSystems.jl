"""
    CSMakie

Module providing Makie.jl plotting functions for ControlSystemsBase.
These functions are loaded when Makie.jl is available through the package extension system.

## Usage
```julia
using ControlSystemsBase, GLMakie  # or CairoMakie
CSMakie.bodeplot(sys)
CSMakie.nyquistplot(sys)
# etc.
```

All functions will throw an informative error if called without Makie.jl loaded.

!!! danger "Experimental"

    The support for plotting with Makie is currently experimental and at any time subject to breaking changes or removal **not** respecting semantic versioning.
"""
module CSMakie

using ..ControlSystemsBase: LTISystem

function bodeplot end
function bodeplot! end

"""
    nyquistplot(sys, [w]; Ms_circles=Float64[], Mt_circles=Float64[], unit_circle=false, critical_point=-1, balance=true, adaptive=true, polar=false, rlimits=(:origin, nothing), kwargs...)

Makie version of [`ControlSystemsBase.nyquistplot`](@ref). In addition to the keyword arguments of the Plots version, this version supports
- `polar`: If `true`, the Nyquist curve is drawn in a `Makie.PolarAxis` with polar grid lines (circles of constant magnitude and rays of constant phase) instead of a Cartesian axis.
- `rlimits`: Radial axis limits, e.g., `(0, 3)`. Only used when `polar = true`. Limiting the radial axis is useful when the system has poles on or close to the imaginary axis, causing the Nyquist curve to reach very large magnitudes.
- `kwargs...`: Any remaining keyword arguments are forwarded to the axis constructor (`Makie.PolarAxis` when `polar = true`, otherwise `Makie.Axis`) and may override the defaults set by this function. This allows customization of, e.g., grid density and appearance:
```julia
CSMakie.nyquistplot(sys; polar=true, rlimits=(0, 3),
    rticks = 0:0.5:3,                                  # magnitude circles
    thetaminorticks = Makie.IntervalsBetween(3),       # phase rays every 15°
    thetaminorgridvisible = true, rminorgridvisible = true)
```
"""
function nyquistplot end
function nyquistplot! end


function sigmaplot end
function sigmaplot! end

function marginplot end
function marginplot! end

function pzmap end
function pzmap! end

function nicholsplot end
function nicholsplot! end

function rgaplot end
function rgaplot! end

function rlocusplot end
function rlocusplot! end

function leadlinkcurve end
function leadlinkcurve! end

# Export all functions and their mutating versions
export bodeplot, bodeplot!, nyquistplot, nyquistplot!, sigmaplot, sigmaplot!, 
       marginplot, marginplot!, pzmap, pzmap!, nicholsplot, nicholsplot!, 
       rgaplot, rgaplot!, rlocusplot, rlocusplot!, leadlinkcurve, leadlinkcurve!

end # module CSMakie