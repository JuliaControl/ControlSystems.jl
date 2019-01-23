struct DelayLtiSystem{T} <: LTISystem
    P::StateSpace{T,Matrix{T}}
    Tau::Vector{Float64} # The length of the vector tau implicitly defines the partitionging of P

    # function DelayLtiSystem(P::StateSpace{T, MT}, Tau::Vector{Float64})
    #     if ControlSystems.noutputs(P) < length(Tau) ||
    #         ControlSystems.noutputs(P) < length(Tau)
    #         error("Length of time-vector is too long given the size of the partitioned system P.")
    #     end
    #     new{T}(P, Tau)
    # end
end

DelayLtiSystem(sys::StateSpace{T,MT}) where {T, MT} = DelayLtiSystem{T}(sys, Float64[])
DelayLtiSystem{T}(sys::StateSpace) where T = DelayLtiSystem{T}(sys, Float64[])


# TODO: Think through these promotions and conversions
Base.promote_rule(::Type{<:StateSpace{T1}}, ::Type{DelayLtiSystem{T2}}) where {T1<:Number,T2<:Number} = DelayLtiSystem{promote_type(T1,T2)}
Base.promote_rule(::Type{AbstractMatrix{T1}}, ::Type{DelayLtiSystem{T2}}) where {T1<:Number,T2<:Number} = DelayLtiSystem{promote_type(T1,T2)}
Base.promote_rule(::Type{T1}, ::Type{DelayLtiSystem{T2}}) where {T1<:Number,T2<:Number} = DelayLtiSystem{promote_type(T1,T2)}
#Base.promote_rule(::Type{<:UniformScaling}, ::Type{S}) where {S<:DelayLtiSystem} = DelayLtiSystem{T}

Base.convert(::Type{DelayLtiSystem{T}}, sys::StateSpace) where T = DelayLtiSystem{T}(sys)
Base.convert(::Type{DelayLtiSystem{T1}}, d::T2) where {T1,T2 <: Number} = DelayLtiSystem{T1}(ss(d))
#Base.convert(::Type{DelayLtiSystem{T}}, sys::DelayLtiSystem) where T = DelayLtiSystem{T}(StateSpace{T,Matrix{T}}(sys))


ninputs(sys::DelayLtiSystem) = ninputs(sys.P) - length(sys.Tau)
noutputs(sys::DelayLtiSystem) = noutputs(sys.P) - length(sys.Tau)
nstates(sys::DelayLtiSystem) = nstates(sys.P)

function +(sys1::DelayLtiSystem, sys2::DelayLtiSystem)
    s1 = PartionedStateSpace(sys1.P, ninputs(sys1), noutputs(sys1))
    s2 = PartionedStateSpace(sys2.P, ninputs(sys2), noutputs(sys2))

    s_new = s1 + s2
    Tau_new = [sys1.Tau; sys2.Tau]

    DelayLtiSystem(s_new.P, Tau_new)
end


function *(sys1::DelayLtiSystem, sys2::DelayLtiSystem)
    s1 = PartionedStateSpace(sys1.P, ninputs(sys1), noutputs(sys1))
    s2 = PartionedStateSpace(sys2.P, ninputs(sys2), noutputs(sys2))

    s_new = s1 * s2
    Tau_new = [sys1.Tau; sys2.Tau]

    DelayLtiSystem(s_new.P, Tau_new)
end


function feedback(sys1::DelayLtiSystem, sys2::DelayLtiSystem)
    s1 = PartionedStateSpace(sys1.P, ninputs(sys1), noutputs(sys1))
    s2 = PartionedStateSpace(sys2.P, ninputs(sys2), noutputs(sys2))

    s_new = feedback(s1, s2)
    Tau_new = [sys1.Tau; sys2.Tau]

    DelayLtiSystem(s_new.P, Tau_new)
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
