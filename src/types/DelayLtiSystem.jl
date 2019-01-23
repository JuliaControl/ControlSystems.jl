struct DelayLtiSystem{T} <: LTISystem
    P::PartionedStateSpace{StateSpace{T,Matrix{T}}}
    Tau::Vector{Float64} # The length of the vector tau implicitly defines the partitionging of P

    # function DelayLtiSystem(P::StateSpace{T, MT}, Tau::Vector{Float64})
    #     if ControlSystems.noutputs(P) < length(Tau) ||
    #         ControlSystems.noutputs(P) < length(Tau)
    #         error("Length of time-vector is too long given the size of the partitioned system P.")
    #     end
    #     new{T}(P, Tau)
    # end
end

# QUESTION: would psys be a good standard variable name for a PartionedStateSpace
#           and perhaps dsys for a delayed system, (ambigous with discrete system though)
function DelayLtiSystem{T}(sys::StateSpace, Tau::Vector{Float64}) where T
    nu = ninputs(sys) - length(Tau)
    ny = noutputs(sys) - length(Tau)

    if nu < 0  || ny < 0
        error("Time vector is too long.")
    end

    psys = PartionedStateSpace(sys, ny, nu)
    DelayLtiSystem{T}(psys, Tau)
end

DelayLtiSystem(sys::StateSpace{T,MT}, Tau::Vector{Float64}) where {T, MT} = DelayLtiSystem{T}(sys, Tau)
DelayLtiSystem(sys::StateSpace{T,MT}) where {T, MT} = DelayLtiSystem{T}(sys, Float64[])
DelayLtiSystem{T}(sys::StateSpace) where T = DelayLtiSystem{T}(sys, Float64[])


# TODO: Think through these promotions and conversions
Base.promote_rule(::Type{<:StateSpace{T1}}, ::Type{DelayLtiSystem{T2}}) where {T1<:Number,T2<:Number} = DelayLtiSystem{promote_type(T1,T2)}
Base.promote_rule(::Type{AbstractMatrix{T1}}, ::Type{DelayLtiSystem{T2}}) where {T1<:Number,T2<:Number} = DelayLtiSystem{promote_type(T1,T2)}
Base.promote_rule(::Type{T1}, ::Type{DelayLtiSystem{T2}}) where {T1<:Number,T2<:Number} = DelayLtiSystem{promote_type(T1,T2)}
#Base.promote_rule(::Type{<:UniformScaling}, ::Type{S}) where {S<:DelayLtiSystem} = DelayLtiSystem{T}

Base.convert(::Type{DelayLtiSystem{T}}, sys::StateSpace) where T = DelayLtiSystem{T}(sys)
Base.convert(::Type{DelayLtiSystem{T1}}, d::T2) where {T1,T2 <: Number} = DelayLtiSystem{T1}(ss(d))
Base.convert(::Type{DelayLtiSystem{T}}, sys::DelayLtiSystem) where T = DelayLtiSystem{T}(StateSpace{T,Matrix{T}}(sys.P.P), sys.Tau)


ninputs(sys::DelayLtiSystem) = size(sys.P.P, 2) - length(sys.Tau)
noutputs(sys::DelayLtiSystem) = size(sys.P.P, 1) - length(sys.Tau)
nstates(sys::DelayLtiSystem) = nstates(sys.P.P)

function +(sys1::DelayLtiSystem, sys2::DelayLtiSystem)
    psys_new = sys1.P + sys2.P
    Tau_new = [sys1.Tau; sys2.Tau]

    DelayLtiSystem(psys_new.P, Tau_new)
end


function *(sys1::DelayLtiSystem, sys2::DelayLtiSystem)
    psys_new = sys1.P * sys2.P
    Tau_new = [sys1.Tau; sys2.Tau]

    DelayLtiSystem(psys_new.P, Tau_new)
end


function feedback(sys1::DelayLtiSystem, sys2::DelayLtiSystem)
    psys_new = feedback(sys1.P, sys2.P)
    Tau_new = [sys1.Tau; sys2.Tau]

    DelayLtiSystem(psys_new.P, Tau_new)
end

function delay(tau::Real, T::Type{<:Number}=Float64)
    return DelayLtiSystem(ControlSystems.ss([zero(T) one(T); one(T) zero(T)]), [float(tau)])
end

# function exp(G::TransferFunction)
#     if (size(G.matrix) != [1, 1]) || ~isone(G.matrix[1].den) || length(G.matrix[1].num) >= 2
#         error("exp only accepts TransferFunction arguemns of the form (a*s + b)")
#     end
#
#     a = G.matrix[1].num[1]
#     b = G.matrix[1].num[0]
#
#     if a > 0
#         error("Delay needs to be causal.")
#     end
#
#     return exp(b) * delay(a)
# end
