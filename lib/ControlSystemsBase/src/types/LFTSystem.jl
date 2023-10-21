abstract type LFTSystem{TE, T} <: LTISystem{TE} end
timeevol(sys::LFTSystem) = timeevol(sys.P)

function *(n::Number, sys::LFTSystem)
    new_C = [sys.P.C1*n; sys.P.C2]
    new_D = [sys.P.D11*n sys.P.D12*n; sys.P.D21 sys.P.D22]
    return basetype(sys)(StateSpace(sys.P.A, sys.P.B, new_C, new_D, sys.P.timeevol), feedback_channel(sys))
end

function *(sys::LFTSystem, n::Number)
    new_B = [sys.P.B1*n sys.P.B2]
    new_D = [sys.P.D11*n sys.P.D12; sys.P.D21*n sys.P.D22]
    return basetype(sys)(StateSpace(sys.P.A, new_B, sys.P.C, new_D, sys.P.timeevol), feedback_channel(sys))
end


function +(sys::LFTSystem{TE,T1}, n::T2) where {TE,T1,T2<:Number}
    T = promote_type(T1,T2)
    if T == T1 # T2 can be stored in sys
        +(sys, T1(n))
    else # We need to upgrade sys
        +(basetype(sys){T}(sys), T(n))
    end
end

# Same numeric type
function +(sys::LFTSystem{<:Any,T}, n::T) where {T<:Number}
    ny, nu = size(sys)
    ssold = sys.P.P
    # Add to direct term from input to output
    new_D = copy(ssold.D)
    new_D[1:ny, 1:nu] .+= n

    pnew = PartitionedStateSpace(StateSpace(ssold.A, ssold.B, ssold.C, new_D, Continuous()), ny, nu)
    basetype(sys)(pnew, feedback_channel(sys))
end

function +(sys1::LFTSystem, sys2::LFTSystem)
    psys_new = sys1.P + sys2.P
    Tau_new = [feedback_channel(sys1); feedback_channel(sys2)]

    promote_type(typeof(sys1), typeof(sys2))(psys_new.P, Tau_new)
end

-(sys1::LFTSystem, sys2::LFTSystem) = +(sys1, -sys2)
-(sys::LFTSystem{<:Any, T}) where {T} = *(sys, T(-1))

# Efficient subtraction with number
-(sys::LFTSystem, n::T) where {T <:Number} = +(sys, -n)
-(n::T, sys::LFTSystem) where {T <:Number} = +(-sys, n)

# Test equality (of realizations)
function ==(sys1::LFTSystem, sys2::LFTSystem)
    all(getfield(sys1, f) == getfield(sys2, f) for f in fieldnames(typeof(sys1)))
end

function *(sys1::LFTSystem{TE, T1}, sys2::LFTSystem{TE, T2}) where {TE, T1, T2}
    psys_new = sys1.P * sys2.P
    Tau_new = [feedback_channel(sys1); feedback_channel(sys2)]
    # T = promote_type(T1, T2)

    basetype(sys1)(psys_new.P, Tau_new)
end


function feedback(sys1::LFTSystem, sys2::LFTSystem; kwargs...)
    isempty(kwargs) || error("The advanced interface to the function `feedback` (with connection keyword arguments) is currently not supported for LFT systems (such as nonlinear and time-delay systems)")
    psys_new = feedback(sys1.P, sys2.P)
    Tau_new = [feedback_channel(sys1); feedback_channel(sys2)]

    promote_type(typeof(sys1), typeof(sys2))(psys_new.P, Tau_new)
end




ninputs(sys::LFTSystem) = size(sys.P.P, 2) - length(feedback_channel(sys))
noutputs(sys::LFTSystem) = size(sys.P.P, 1) - length(feedback_channel(sys))
nstates(sys::LFTSystem) = nstates(sys.P.P)

Base.size(sys::LFTSystem) = (noutputs(sys), ninputs(sys))

# Fallbacks, TODO We should sort this out for all types, maybe after SISO/MIMO
# {Array, Number}, Colon
Base.getindex(sys::LFTSystem, i, ::Colon) =
    getindex(sys, index2range(i), 1:size(sys,2))
# Colon, {Array, Number}
Base.getindex(sys::LFTSystem, ::Colon, j) =
    getindex(sys, 1:size(sys,1), index2range(j))
Base.getindex(sys::LFTSystem, ::Colon, ::Colon) =
    getindex(sys, 1:size(sys,1), 1:size(sys,2))
# Should just catch Number, Number, or Colon, Colon
Base.getindex(sys::LFTSystem, i, j) =
    getindex(sys, index2range(i), index2range(j))

function Base.getindex(sys::LFTSystem, i::AbstractArray, j::AbstractArray)
    ny, nu = size(sys)
    # Can't use "boundscheck" since not AbstractArray
    imin, imax = extrema(i)
    jmin, jmax = extrema(j)
    if imax > ny || imin < 1 || jmax > nu || jmin < 1
        throw(BoundsError(sys, (i,j)))
    end
    nrow, ncol = size(sys.P.P)
    rowidx = [i; (ny+1):nrow] # Output plus feedback terms
    colidx = [j; (nu+1):ncol] # Input plus feedback terms
    typeof(sys)(StateSpace(
        sys.P.A[:,      :],
        sys.P.B[:,      colidx],
        sys.P.C[rowidx, :],
        sys.P.D[rowidx, colidx],
        sys.P.timeevol), feedback_channel(sys))
end

function append(systems::LFTT...) where LFTT <: LFTSystem
    timeevol = common_timeevol(systems...)
    A   = blockdiag(s.P.A for s in systems)
    B1  = blockdiag(s.P.B1 for s in systems)
    B2  = blockdiag(s.P.B2 for s in systems)
    C1  = blockdiag(s.P.C1 for s in systems)
    C2  = blockdiag(s.P.C2 for s in systems)
    D11 = blockdiag(s.P.D11 for s in systems)
    D12 = blockdiag(s.P.D12 for s in systems)
    D21 = blockdiag(s.P.D21 for s in systems)
    D22 = blockdiag(s.P.D22 for s in systems)
    LFTT(
        PartitionedStateSpace(A,B1,B2,C1,C2,D11,D12,D21,D22,timeevol),
        reduce(vcat, [feedback_channel(s) for s in systems])
    )
end

