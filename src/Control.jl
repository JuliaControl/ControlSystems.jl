module Control

export  LTISystem,
        StateSpace,
        TransferFunction,
        ss,
        tf,
        # Linear Algebra
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
        # Model Simplification
        sminreal,
        # Stability Analysis
        pole,
        tzero,
        gain,
        zpkdata,
        damp,
        dampreport,
        markovparam,
        # Connections
        append,
        series,
        parallel,
        # Discrete
        c2d,
        # Time Response
        step,
        impulse,
        lsim,
        # Frequency Response
        freqresp,
        evalfr,
        bode,
        nyquist,
        sigma

using Requires

import Base: +, -, *, /, (./), (==), (.+), (.-), (.*)

include("types/lti.jl")
include("types/transferfunction.jl")
include("types/statespace.jl")
include("types/tf2ss.jl")

include("connections.jl")
include("discrete.jl")
include("matrix_comps.jl")
include("simplification.jl")
include("synthesis.jl")
include("analysis.jl")
include("timeresp.jl")
include("freqresp.jl")
include("utilities.jl")

# The path has to be evaluated upon initial import
const __CONTROL_SOURCE_DIR__ = dirname(Base.source_path())
@require PyPlot include(joinpath(__CONTROL_SOURCE_DIR__, "plotting.jl"))

end
