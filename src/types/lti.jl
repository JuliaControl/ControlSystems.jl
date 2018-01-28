abstract type LTISystem end
+(sys1::LTISystem, sys2::LTISystem) = +(promote(sys1, sys2)...)
-(sys1::LTISystem, sys2::LTISystem) = -(promote(sys1, sys2)...)
*(sys1::LTISystem, sys2::LTISystem) = *(promote(sys1, sys2)...)
/(sys1::LTISystem, sys2::LTISystem) = /(promote(sys1, sys2)...)

const AbstractNumberVector = AbstractVector{<:Number}
const AbstractNumberMatrix = AbstractMatrix{<:Number}
const BlasNumber = Union{Base.LinAlg.BlasFloat, Base.LinAlg.BlasInt}
