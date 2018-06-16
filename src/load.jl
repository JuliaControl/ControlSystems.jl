using Polynomials

const BlasNumber = Union{Base.LinAlg.BlasFloat, Base.LinAlg.BlasInt}

import Base: +, -, *, /, (./), (==), (.+), (.-), (.*), (!=), isapprox, convert, promote_op, num, den


#isapprox(x::Poly, y::Poly; kwargs...) = isapprox(promote(x, y)...; kwargs...)

abstract type AbstractSystem end

include("types/Lti.jl")

abstract type SisoTf{T<:Number} end


# Transfer functions and tranfer function elemements
include("types/TransferFunction.jl")
include("types/SisoTfTypes/SisoZpk.jl")
include("types/SisoTfTypes/SisoRational.jl")
include("types/SisoTfTypes/SisoGeneralized.jl")
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

using OrdinaryDiffEq
include("simulators.jl")
#include("pid_design.jl")


using Base.Test
import Base.isapprox
include("../test/framework.jl")


#function Base.blkdiag(systems::Union{LTISystem,Number,AbstractMatrix{<:Number}}...)
#end
