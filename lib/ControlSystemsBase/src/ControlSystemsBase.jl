module ControlSystemsBase

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
        isproper,
        StaticStateSpace,
        to_static,
        to_sized,
        system_name,
        input_names,
        output_names,
        state_names,
        # Linear Algebra
        balance,
        balance_statespace,
        are,
        lqr,
        kalman,
        covar,
        norm,
        hinfnorm,
        linfnorm,
        gram,
        grampd,
        ctrb,
        obsv,
        controllability,
        observability,
        place,
        place_knvd,
        # Model Simplification
        reduce_sys,
        sminreal,
        minreal,
        balreal,
        baltrunc,
        similarity_transform,
        time_scale,
        innovation_form,
        observer_predictor,
        observer_filter,
        observer_controller,
        # Stability Analysis
        isstable,
        poles,
        tzeros,
        dcgain,
        zpkdata,
        damp,
        dampreport,
        markovparam,
        margin,
        delaymargin,
        gangoffour,
        extended_gangoffour,
        relative_gain_array,
        # Connections
        append,
        series,
        parallel,
        array2mimo,
        feedback,
        feedback2dof,
        starprod,
        lft,
        sensitivity,
        input_sensitivity,
        output_sensitivity,
        comp_sensitivity,
        input_comp_sensitivity,
        output_comp_sensitivity,
        G_PS,
        G_CS,
        # Discrete
        c2d,
        c2d_x0map,
        d2c,
        # Time Response
        step,
        impulse,
        lsim,
        lsim!,
        LsimWorkspace,
        stepinfo,
        StepInfo,
        # Frequency Response
        freqresp, freqrespv, freqresp!,
        evalfr,
        bode, bodev,
        bodemag!,
        BodemagWorkspace,
        nyquist, nyquistv,
        sigma, sigmav,
        # delay systems
        delay,
        pade,
        nonlinearity,
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
        ssdata,
        add_input,
        add_output


# QUESTION: are these used? LaTeXStrings, Requires, IterTools
using RecipesBase, LaTeXStrings, LinearAlgebra
import Polynomials
import Polynomials: Polynomial, coeffs
import Base: +, -, *, /, (==), (!=), isapprox, convert, promote_op
import Base: getproperty, getindex
import Base: exp # for exp(-s)
import LinearAlgebra: BlasFloat

export lyap # Make sure LinearAlgebra.lyap is available
import Printf
import Printf: @printf, @sprintf
import DSP
import DSP: conv
using ForwardDiff
import MatrixPencils
using MacroTools
using MatrixEquations
using UUIDs # to load Plots in gangoffourplot
using StaticArrays, Polyester

abstract type AbstractSystem end


include("types/result_types.jl")
include("types/TimeEvolution.jl")
## Added interface:
#   timeevol(Lti) -> TimeEvolution (not exported)


include("types/Lti.jl")

include("types/SisoTf.jl")

# Transfer functions and transfer function elements
include("types/TransferFunction.jl")
include("types/SisoTfTypes/SisoZpk.jl")
include("types/SisoTfTypes/SisoRational.jl")
include("types/SisoTfTypes/promotion.jl")
include("types/SisoTfTypes/conversion.jl")

include("types/StateSpace.jl")

# TODO Sample time
include("types/PartitionedStateSpace.jl")
include("types/LFTSystem.jl")
include("types/DelayLtiSystem.jl")
include("types/HammersteinWiener.jl")

# Convenience constructors
include("types/tf.jl")
include("types/zpk.jl")

include("utilities.jl")

include("types/promotion.jl")
include("types/conversion.jl")
include("connections.jl")
include("sensitivity_functions.jl")

# Analysis
include("freqresp.jl")
include("timeresp.jl")

include("matrix_comps.jl")
include("simplification.jl")

include("discrete.jl")
include("analysis.jl")
include("synthesis.jl")

include("pid_design.jl")

include("demo_systems.jl")

include("delay_systems.jl")
include("hammerstein_weiner.jl")
include("nonlinear_components.jl")

include("types/staticsystems.jl")

include("plotting.jl")
include("dsp.jl")

@deprecate pole poles
@deprecate tzero tzeros
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
const __ControlSystemsBase_SOURCE_DIR__ = dirname(Base.source_path())

function __init__()
    Base.Experimental.register_error_hint(MethodError) do io, exc, argtypes, kwargs
        if exc.f ∈ (lsim, step, impulse) && argtypes[1] <: LTISystem{Continuous}
            print(io, "\n$(exc.f) with continuous-time systems, including delay systems and nonlinear systems, require the user to first ")
            printstyled(io, "install and load ControlSystems.jl, or pass the keyword method = :zoh", color=:green, bold=true)
            print(io, " for automatic discretization (applicable to systems without delays or nonlinearities only).")
        end
        plots_id = Base.PkgId(UUID("91a5bcdd-55d7-5caf-9e0b-520d859cae80"), "Plots")
        if nameof(exc.f) === :plot && parentmodule(argtypes[1]) == @__MODULE__() && !haskey(Base.loaded_modules, plots_id)
            printstyled(io, "\nPlotting is not available unless Plots.jl is loaded manually. Call `using Plots` before plotting.", color=:green, bold=true)
        elseif (exc.f == /) && argtypes[2] <: DelayLtiSystem
            print(io, "A delayed system can not be inverted. Consider use of the function `feedback`.")
        end
    end
end

include("precompilation.jl")

end
