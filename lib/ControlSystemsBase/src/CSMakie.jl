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