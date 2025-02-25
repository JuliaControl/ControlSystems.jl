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
        stab_unstab,
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
        resolvent,
        input_resolvent,
        # Discrete
        c2d,
        c2d_x0map,
        d2c,
        d2c_exact,
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
        thiran,
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


using RecipesBase, LinearAlgebra
import Polynomials
import Polynomials: Polynomial, coeffs
import Base: +, -, *, /, (==), (!=), isapprox, convert, promote_op
import Base: getproperty, getindex
import Base: exp # for exp(-s)
import LinearAlgebra: BlasFloat

export lyap # Make sure LinearAlgebra.lyap is available
export plyap
import Printf
import Printf: @printf, @sprintf
import Polynomials: conv # TODO: replace this internal function with something public
using ForwardDiff
import MatrixPencils
using MacroTools
using MatrixEquations
using UUIDs # to load Plots in gangoffourplot
using StaticArraysCore

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

@deprecate pole poles
@deprecate tzero tzeros
@deprecate num numvec
@deprecate den denvec
@deprecate norminf hinfnorm
@deprecate diagonalize(s::AbstractStateSpace, digits) diagonalize(s::AbstractStateSpace)
@deprecate luenberger(sys, p) place(sys, p, :o)
@deprecate luenberger(A, C, p) place(A, C, p, :o)
# There are some deprecations in pid_control.jl for laglink/leadlink/leadlinkat

"""
    Gs, k = seriesform(G::TransferFunction{Discrete})

Convert a transfer function `G` to a vector of second-order transfer functions and a scalar gain `k`, the product of which equals `G`.

!!! note
    This function requires the user to load the package DSP.jl.
"""
seriesform(a) = error(a isa TransferFunction{<:Discrete} ? "seriesform requires the user to load the package DSP" : "seriesform requires a discrete-time TransferFunction (and the package DSP.jl to be loaded)")

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
        elseif exc.f ∈ (eigvals!, ) && argtypes[1] <: AbstractMatrix{<:Number}
            printstyled(io, "\nComputing eigenvalues of a matrix with exotic element types may require `using GenericSchur`.", color=:green, bold=true)
        end
        plots_id = Base.PkgId(UUID("91a5bcdd-55d7-5caf-9e0b-520d859cae80"), "Plots")
        if exc.f isa Function && nameof(exc.f) === :plot && parentmodule(argtypes[1]) == @__MODULE__() && !haskey(Base.loaded_modules, plots_id)
            printstyled(io, "\nPlotting is not available unless Plots.jl is loaded manually. Call `using Plots` before plotting.", color=:green, bold=true)
        elseif (exc.f == /) && argtypes[2] <: DelayLtiSystem
            print(io, "A delayed system can not be inverted. Consider use of the function `feedback`.")
        end
    end
end

include("precompilation.jl")

end
