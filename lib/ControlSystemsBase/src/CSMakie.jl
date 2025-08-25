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
"""
module CSMakie

using ..ControlSystemsBase: LTISystem

# Error message for when Makie is not loaded
const MAKIE_NOT_LOADED_ERROR = """
Makie.jl must be loaded to use CSMakie plotting functions.
Please run one of:
  using GLMakie      # for interactive plots
  using CairoMakie   # for static plots
  using WGLMakie     # for web-based plots
"""

# Define stub functions that will be overloaded by the extension

"""
    bodeplot(sys::LTISystem, w=nothing; kwargs...)
    bodeplot(sys::Vector{<:LTISystem}, w=nothing; kwargs...)

Create a Bode plot using Makie.jl. Requires Makie to be loaded.

See the main `bodeplot` documentation for available options.
"""
function bodeplot(args...; kwargs...)
    error(MAKIE_NOT_LOADED_ERROR)
end

"""
    nyquistplot(sys::LTISystem, w=nothing; kwargs...)
    nyquistplot(sys::Vector{<:LTISystem}, w=nothing; kwargs...)

Create a Nyquist plot using Makie.jl. Requires Makie to be loaded.

See the main `nyquistplot` documentation for available options.
"""
function nyquistplot(args...; kwargs...)
    error(MAKIE_NOT_LOADED_ERROR)
end

"""
    sigmaplot(sys::LTISystem, w=nothing; kwargs...)
    sigmaplot(sys::Vector{<:LTISystem}, w=nothing; kwargs...)

Create a singular value plot using Makie.jl. Requires Makie to be loaded.

See the main `sigmaplot` documentation for available options.
"""
function sigmaplot(args...; kwargs...)
    error(MAKIE_NOT_LOADED_ERROR)
end

"""
    marginplot(sys::LTISystem, w=nothing; kwargs...)
    marginplot(sys::Vector{<:LTISystem}, w=nothing; kwargs...)

Create a margin plot using Makie.jl. Requires Makie to be loaded.

See the main `marginplot` documentation for available options.
"""
function marginplot(args...; kwargs...)
    error(MAKIE_NOT_LOADED_ERROR)
end

"""
    pzmap(sys::LTISystem; kwargs...)
    pzmap(sys::Vector{<:LTISystem}; kwargs...)

Create a pole-zero map using Makie.jl. Requires Makie to be loaded.

See the main `pzmap` documentation for available options.
"""
function pzmap(args...; kwargs...)
    error(MAKIE_NOT_LOADED_ERROR)
end

"""
    pzmap!(ax, sys::LTISystem; kwargs...)
    pzmap!(ax, sys::Vector{<:LTISystem}; kwargs...)

Add a pole-zero map to an existing Makie axis. Requires Makie to be loaded.
"""
function pzmap!(args...; kwargs...)
    error(MAKIE_NOT_LOADED_ERROR)
end

"""
    nicholsplot(sys::LTISystem, w=nothing; kwargs...)
    nicholsplot(sys::Vector{<:LTISystem}, w=nothing; kwargs...)

Create a Nichols chart using Makie.jl. Requires Makie to be loaded.

See the main `nicholsplot` documentation for available options.
"""
function nicholsplot(args...; kwargs...)
    error(MAKIE_NOT_LOADED_ERROR)
end

"""
    rgaplot(sys::LTISystem, w=nothing; kwargs...)
    rgaplot(sys::Vector{<:LTISystem}, w=nothing; kwargs...)

Create a relative gain array plot using Makie.jl. Requires Makie to be loaded.

See the main `rgaplot` documentation for available options.
"""
function rgaplot(args...; kwargs...)
    error(MAKIE_NOT_LOADED_ERROR)
end

"""
    rlocusplot(P::LTISystem, K=500; kwargs...)

Create a root locus plot using Makie.jl. Requires Makie to be loaded.

See the main `rlocusplot` documentation for available options.
"""
function rlocusplot(args...; kwargs...)
    error(MAKIE_NOT_LOADED_ERROR)
end

"""
    leadlinkcurve(start=1; kwargs...)

Plot the phase advance curve for a lead link using Makie.jl. Requires Makie to be loaded.

See the main `leadlinkcurve` documentation for available options.
"""
function leadlinkcurve(args...; kwargs...)
    error(MAKIE_NOT_LOADED_ERROR)
end

"""
    setPlotScale(str)

Set the default scale of magnitude in `bodeplot` and `sigmaplot` for Makie plots.
`str` should be either `"dB"` or `"log10"`. The default scale is `"log10"`.

Requires Makie to be loaded.
"""
function setPlotScale(str::AbstractString)
    error(MAKIE_NOT_LOADED_ERROR)
end

# Also export the mutating versions that will be defined by the extension
export bodeplot, nyquistplot, sigmaplot, marginplot, pzmap, pzmap!,
       nicholsplot, rgaplot, rlocusplot, leadlinkcurve, setPlotScale

end # module CSMakie