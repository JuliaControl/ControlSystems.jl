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
function DelayLtiSystem{T}(sys::StateSpace, Tau::Vector{Float64}) where T<:Number
    nu = ninputs(sys) - length(Tau)
    ny = noutputs(sys) - length(Tau)

    if nu < 0  || ny < 0
        throw(ArgumentError("Time vector is too long."))
    end

    psys = PartionedStateSpace{StateSpace{T,Matrix{T}}}(sys, nu, ny)
    DelayLtiSystem{T}(psys, Tau)
end
# For converting DelayLtiSystem{T} to different T
DelayLtiSystem{T}(sys::DelayLtiSystem) where T = DelayLtiSystem{T}(PartionedStateSpace{StateSpace{T,Matrix{T}}}(sys.P), Float64[])
DelayLtiSystem{T}(sys::StateSpace) where T = DelayLtiSystem{T}(sys, Float64[])

# From StateSpace, infer type
DelayLtiSystem(sys::StateSpace{T,MT}, Tau::Vector{Float64}) where {T, MT} = DelayLtiSystem{T}(sys, Tau)
DelayLtiSystem(sys::StateSpace{T,MT}) where {T, MT} = DelayLtiSystem{T}(sys, Float64[])

# From TransferFunction, infer type TODO Use proper constructor instead of convert here when defined
DelayLtiSystem(sys::TransferFunction{S}) where {T,S<:SisoTf{T}} = DelayLtiSystem(convert(StateSpace{T, Matrix{T}}, sys))

# TODO: Think through these promotions and conversions
Base.promote_rule(::Type{AbstractMatrix{T1}}, ::Type{DelayLtiSystem{T2}}) where {T1<:Number,T2<:Number} = DelayLtiSystem{promote_type(T1,T2)}
Base.promote_rule(::Type{T1}, ::Type{DelayLtiSystem{T2}}) where {T1<:Number,T2<:Number} = DelayLtiSystem{promote_type(T1,T2)}
#Base.promote_rule(::Type{<:UniformScaling}, ::Type{S}) where {S<:DelayLtiSystem} = DelayLtiSystem{T}

Base.convert(::Type{DelayLtiSystem{T}}, sys::StateSpace) where T = DelayLtiSystem{T}(sys)
Base.convert(::Type{DelayLtiSystem{T1}}, d::T2) where {T1,T2 <: Number} = DelayLtiSystem{T1}(ss(d))
# Catch convertsion between T
Base.convert(::Type{S}, sys::DelayLtiSystem)  where {T, S<:DelayLtiSystem{T}} =
    sys isa S ? sys : S(StateSpace{T,Matrix{T}}(sys.P.P), sys.Tau)

Base.convert(::Type{DelayLtiSystem{T}}, sys::TransferFunction) where T = DelayLtiSystem{T}(convert(StateSpace{T, Matrix{T}}, sys))


function *(sys::DelayLtiSystem, n::Number)
    new_C = [sys.P.C1*n; sys.P.C2]
    new_D = [sys.P.D11*n sys.P.D12; sys.P.D21*n sys.P.D22]
    return DelayLtiSystem(StateSpace(sys.P.A, sys.P.B, new_C, new_D, sys.P.Ts), sys.Tau)
end
*(n::Number, sys::DelayLtiSystem) = *(sys, n)

function +(sys::DelayLtiSystem{T1}, n::T2) where {T1,T2<:Number}
    T = promote_type(T1,T2)
    if T == T1 # T2 can be stored in sys
        +(sys, T1(n))
    else # We need to upgrade sys
        +(DelayLtiSystem{T}(sys), T(n))
    end
end

# Efficient addition
function +(sys::DelayLtiSystem{T}, n::T) where {T<:Number}
    ny, nu = size(sys)
    ssold = sys.P.P
    # Add to direct term from input to output
    new_D = ssold.D
    new_D[1:ny, 1:nu] .+= n

    pnew = PartionedStateSpace(StateSpace(ssold.A, ssold.B, ssold.C, new_D, 0.0), ny, nu)
    DelayLtiSystem(pnew, sys.Tau)
end

# Efficient subtraction with number
-(sys::DelayLtiSystem, n::T) where {T <:Number} = +(sys, -n)
-(n::T, sys::DelayLtiSystem) where {T <:Number} = +(-sys, n)

ninputs(sys::DelayLtiSystem) = size(sys.P.P, 2) - length(sys.Tau)
noutputs(sys::DelayLtiSystem) = size(sys.P.P, 1) - length(sys.Tau)
nstates(sys::DelayLtiSystem) = nstates(sys.P.P)

Base.size(sys::DelayLtiSystem) = (noutputs(sys), ninputs(sys))

# Fallbacks, TODO We should sort this out for all types, maybe after SISO/MIMO
# {Array, Number}, Colon
Base.getindex(sys::DelayLtiSystem, i, ::Colon) =
    getindex(sys, index2range(i), 1:size(sys,2))
# Colon, {Array, Number}
Base.getindex(sys::DelayLtiSystem, ::Colon, j) =
    getindex(sys, 1:size(sys,1), index2range(j))
Base.getindex(sys::DelayLtiSystem, ::Colon, ::Colon) =
    getindex(sys, 1:size(sys,1), 1:size(sys,2))
# Should just catch Number, Number, or Colon, Colon
Base.getindex(sys::DelayLtiSystem, i, j) =
    getindex(sys, index2range(i), index2range(j))

function Base.getindex(sys::DelayLtiSystem, i::AbstractArray, j::AbstractArray)
    ny, nu = size(sys)
    # Cant use "boundscheck" since not AbstractArray
    imin, imax = extrema(i)
    jmin, jmax = extrema(j)
    if imax > ny || imin < 1 || jmax > nu || jmin < 1
        throw(BoundsError(sys, (i,j)))
    end
    nrow, ncol = size(sys.P.P)
    rowidx = [collect(i); collect((ny+1):nrow)] # Output plus delay terms
    colidx = [collect(j); collect((nu+1):ncol)] # Input plus delay terms
    DelayLtiSystem(StateSpace(
        sys.P.A[:,      :],
        sys.P.B[:,      colidx],
        sys.P.C[rowidx, :],
        sys.P.D[rowidx, colidx],
        sys.P.Ts), sys.Tau)
end

function +(sys1::DelayLtiSystem, sys2::DelayLtiSystem)
    psys_new = sys1.P + sys2.P
    Tau_new = [sys1.Tau; sys2.Tau]

    DelayLtiSystem(psys_new.P, Tau_new)
end

-(sys1::DelayLtiSystem, sys2::DelayLtiSystem) = +(sys1, -sys2)
-(sys::DelayLtiSystem{T}) where {T} = *(sys, T(-1))


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
