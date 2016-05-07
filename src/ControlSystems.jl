module ControlSystems

export  LTISystem,
        StateSpace,
        TransferFunction,
        ss,
        tf,
        tfg,
        zpk,
        ss2tf,
        # Linear Algebra
        balance,
        care,
        dare,
        dlyap,
        lqr,
        dlqr,
        kalman,
        dkalman,
        covar,
        norm,
        gram,
        ctrb,
        obsv,
        place,
        # Model Simplification
        reduce_sys,
        sminreal,
        minreal,
        balreal,
        # Stability Analysis
        isstable,
        pole,
        tzero,
        dcgain,
        zpkdata,
        damp,
        dampreport,
        markovparam,
        margin,
        delaymargin,
        gangoffour,
        # Connections
        append,
        series,
        parallel,
        feedback,
        feedback2dof,
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
        sigma,
        # utilities
        numpoly,
        denpoly,
        numvec,
        denvec

import Plots
import Base: +, -, *, /, (./), (==), (.+), (.-), (.*), (!=), isapprox, call, convert

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
const __CONTROLSYSTEMS_SOURCE_DIR__ = dirname(Base.source_path())

end
