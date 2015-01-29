module Control

export LTISystem,
       StateSpace,
       TransferFunction,
       ss,
       tf,
       care,
       dare,
       dlyap,
       lqr,
       dlqr,
       covar,
       norm,
       gram,
       ctrb,
       obsv,
       pole,
       tzero,
       damp,
       dampreport,
       append,
       series,
       parallel,
       c2d,
       step,
       impulse,
       lsim

using Requires

include("types/lti.jl")
include("types/transferfunction.jl")
include("types/statespace.jl")
include("types/tf2ss.jl")

include("connections.jl")
include("discrete.jl")
include("matrix_comps.jl")
include("synthesis.jl")
include("analysis.jl")
include("timeresp.jl")
include("utilities.jl")

# The path has to be evaluated upon initial import
const __CONTROL_SOURCE_DIR__ = dirname(Base.source_path())
@require PyPlot include(joinpath(__CONTROL_SOURCE_DIR__, "plotting.jl"))

end
