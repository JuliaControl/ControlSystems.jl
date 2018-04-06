using Polynomials

const BlasNumber = Union{Base.LinAlg.BlasFloat, Base.LinAlg.BlasInt}

import Base: +, -, *, /, (./), (==), (.+), (.-), (.*), (!=), isapprox, convert, promote_op, num, den

include("types/Lti.jl")

abstract type SisoTf{T<:Number} end



include("types/TransferFunction.jl")
include("types/SisoTfTypes/SisoZpk.jl")
include("types/SisoTfTypes/SisoRational.jl")
include("types/SisoTfTypes/SisoGeneralized.jl")
include("types/SisoTfTypes/promotion.jl")
include("types/SisoTfTypes/conversion.jl")

include("types/tf.jl")
include("types/zpk.jl")

include("types/StateSpace.jl")

include("utilities.jl")


include("types/promotion.jl")
include("types/conversion.jl")
include("connections.jl")

# Analysis
include("freqresp.jl")


using Base.Test
import Base.isapprox
include("../../test/framework.jl")
