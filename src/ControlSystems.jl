module ControlSystems

export  LTISystem,
        AbstractStateSpace,
        StateSpace,
        HeteroStateSpace,
        TransferFunction,
        DelayLtiSystem,
        ss,
        tf,
        zpk,
        ss2tf,
        LQG,
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
        similarity_transform,
        innovation_form,
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
        # delay systems
        delay,
        # utilities
        num,    #Deprecated
        den,    #Deprecated
        numvec,
        denvec,
        numpoly,
        denpoly


# QUESTION: are these used? LaTeXStrings, Requires, IterTools
using Polynomials, Plots, LaTeXStrings, LinearAlgebra
using OrdinaryDiffEq, DelayDiffEq
export Plots
import Base: +, -, *, /, (==), (!=), isapprox, convert, promote_op
import Base: getproperty
import LinearAlgebra: BlasFloat
export lyap # Make sure LinearAlgebra.lyap is available
import Printf, Colors
import DSP: conv

abstract type AbstractSystem end

include("types/Lti.jl")

include("types/SisoTf.jl")

# Transfer functions and tranfer function elemements
include("types/TransferFunction.jl")
include("types/SisoTfTypes/polyprint.jl")
include("types/SisoTfTypes/SisoZpk.jl")
include("types/SisoTfTypes/SisoRational.jl")
include("types/SisoTfTypes/promotion.jl")
include("types/SisoTfTypes/conversion.jl")

include("types/StateSpace.jl")

include("types/PartionedStateSpace.jl")
include("types/DelayLtiSystem.jl")

# Convenience constructors
include("types/tf.jl")
include("types/zpk.jl")

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

include("delay_systems.jl")

include("plotting.jl")

@deprecate num numvec
@deprecate den denvec

# The path has to be evaluated upon initial import
const __CONTROLSYSTEMS_SOURCE_DIR__ = dirname(Base.source_path())

end
