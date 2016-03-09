module Control

export  LTISystem,
        StateSpace,
        TransferFunction,
        ss,
        tf,
        zpk,
        ss2tf,
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
        place,
        # Model Simplification
        sminreal,
        minreal,
        # Stability Analysis
        pole,
        tzero,
        gain,
        dcgain,
        zpkdata,
        damp,
        dampreport,
        markovparam,
        margin,
        gangoffour,
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

import Plots

import Base: +, -, *, /, (./), (==), (.+), (.-), (.*), call

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
include("plotting.jl")
include("pid_design.jl")

# The path has to be evaluated upon initial import
const __CONTROL_SOURCE_DIR__ = dirname(Base.source_path())

end
