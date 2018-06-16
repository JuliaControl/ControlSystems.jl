module ControlSystems

const BlasNumber = Union{Base.LinAlg.BlasFloat, Base.LinAlg.BlasInt} # QUESTION: WHy include BlasInt, why not include BlasComplex?

export  LTISystem,
        StateSpace,
        TransferFunction,
        ss,
        tf,
        zpk,
        ss2tf,
        LQG,
        primitivetype,
        # Linear Algebra
        balance,
        care,
        dare,
        dlyap,
        lqr,
        dlqr,
        kalman,
        dkalman,
        lqg,
        lqgi,
        covar,
        norm,
        norminf,
        gram,
        ctrb,
        obsv,
        place,
        # Model Simplification
        reduce_sys,
        sminreal,
        minreal,
        balreal,
        baltrunc,
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
        solve,
        Simulator,
        # Frequency Response
        freqresp,
        evalfr,
        bode,
        nyquist,
        sigma,
        # utilities
        numpoly,
        denpoly


# QUESTION: are these used? LaTeXStrings, Requires, IterTools
using Polynomials, OrdinaryDiffEq, Plots
import Base: +, -, *, /, (./), (==), (.+), (.-), (.*), (!=), isapprox, convert, promote_op

abstract type AbstractSystem end

include("types/Lti.jl")

abstract type SisoTf{T<:Number} end

# Transfer functions and tranfer function elemements
include("types/TransferFunction.jl")
include("types/SisoTfTypes/SisoZpk.jl")
include("types/SisoTfTypes/SisoRational.jl")
include("types/SisoTfTypes/promotion.jl")
include("types/SisoTfTypes/conversion.jl")

include("types/StateSpace.jl")

# Convenience constructors
include("types/tf.jl")
include("types/zpk.jl")
include("types/ss.jl")

include("types/lqg.jl") # QUESTION: is it really motivated to have an LQG type?

include("utilities.jl")

include("types/promotion.jl")
include("types/conversion.jl")
include("connections.jl")

# Analysis
include("freqresp.jl")
include("timeresp.jl")

include("matrix_comps.jl")
include("simplification.jl")

include("discrete.jl")
include("analysis.jl")
include("synthesis.jl")

include("simulators.jl")
include("pid_design.jl")



# The path has to be evaluated upon initial import
const __CONTROLSYSTEMS_SOURCE_DIR__ = dirname(Base.source_path())

end
