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

DelayLtiSystem(sys::StateSpace) = DelayLtiSystem(sys, Float64[])
#Base.convert(::Type{S}, ) where {S<:DelayLtiSystem} =  # Need to handle numerical arguments

Base.promote_rule(::Type{<:StateSpace}, ::Type{S}) where {S<:DelayLtiSystem} = S

ninputs(sys::DelayLtiSystem) = size(sys.P)[2] - length(sys.Tau)
noutputs(sys::DelayLtiSystem) = size(sys.P)[1] - length(sys.Tau)


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
