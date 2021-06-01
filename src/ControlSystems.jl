module ControlSystems

export  LTISystem,
        AbstractStateSpace,
        StateSpace,
        HeteroStateSpace,
        TransferFunction,
        DelayLtiSystem,
        Continuous,
        Discrete,
        ss,
        tf,
        zpk,
        LQG,
        isproper,
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
        hinfnorm,
        linfnorm,
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
        prescale,
        innovation_form,
        observer_predictor,
        observer_controller,
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
        starprod,
        lft,
        # Discrete
        c2d,
        c2d_x0map,
        d2c,
        # Time Response
        step,
        impulse,
        lsim,
        solve,
        Simulator,
        # Frequency Response
        freqresp, freqrespv,
        evalfr,
        bode, bodev,
        nyquist, nyquistv,
        sigma, sigmav,
        # delay systems
        delay,
        pade,
        # demo systems
        ssrand,
        DemoSystems, # A module containing some example systems
        # utilities
        num,    #Deprecated
        den,    #Deprecated
        numvec,
        denvec,
        numpoly,
        denpoly,
        iscontinuous,
        isdiscrete,
        ssdata


# QUESTION: are these used? LaTeXStrings, Requires, IterTools
using Plots, LaTeXStrings, LinearAlgebra
import Polynomials
import Polynomials: Polynomial, coeffs
using OrdinaryDiffEq
export Plots
import Base: +, -, *, /, (==), (!=), isapprox, convert, promote_op
import Base: getproperty, getindex
import Base: exp # for exp(-s)
import LinearAlgebra: BlasFloat
export lyap # Make sure LinearAlgebra.lyap is available
import Printf, Colors
import DSP: conv
import DiffEqCallbacks: SavingCallback, SavedValues
using DelayDiffEq
using MacroTools

abstract type AbstractSystem end

include("types/TimeEvolution.jl")
## Added interface:
#   timeevol(Lti) -> TimeEvolution (not exported)


include("types/Lti.jl")

include("types/SisoTf.jl")

# Transfer functions and tranfer function elemements
include("types/TransferFunction.jl")
include("types/SisoTfTypes/SisoZpk.jl")
include("types/SisoTfTypes/SisoRational.jl")
include("types/SisoTfTypes/promotion.jl")
include("types/SisoTfTypes/conversion.jl")

include("types/StateSpace.jl")

# TODO Sample time
include("types/PartionedStateSpace.jl")
include("types/DelayLtiSystem.jl")

# Convenience constructors
include("types/tf.jl")
include("types/zpk.jl")

include("types/lqg.jl") # QUESTION: is it really motivated to have an LQG type?

include("utilities.jl")
include("auto_grid.jl")

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

include("demo_systems.jl")

include("delay_systems.jl")

include("plotting.jl")

@deprecate num numvec
@deprecate den denvec
@deprecate norminf hinfnorm
@deprecate diagonalize(s::AbstractStateSpace, digits) diagonalize(s::AbstractStateSpace)
@deprecate luenberger(sys, p) place(sys, p, :o)
@deprecate luenberger(A, C, p) place(A, C, p, :o)
# There are some deprecations in pid_control.jl for laglink/leadlink/leadlinkat

function covar(D::Union{AbstractMatrix,UniformScaling}, R)
    @warn "This call is deprecated due to ambiguity, use covar(ss(D), R) or covar(ss(D, Ts), R) instead"
    D*R*D'
end

# The path has to be evaluated upon initial import
const __CONTROLSYSTEMS_SOURCE_DIR__ = dirname(Base.source_path())

end
