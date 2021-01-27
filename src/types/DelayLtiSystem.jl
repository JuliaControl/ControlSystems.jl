"""
    struct DelayLtiSystem{T, S <: Real} <: LTISystem

Represents an LTISystem with internal time-delay. See `?delay` for a convenience constructor.
"""
struct DelayLtiSystem{T,S<:Real} <: LTISystem
    P::PartionedStateSpace{StateSpace{Continuous,T}}
    Tau::Vector{S} # The length of the vector tau implicitly defines the partitionging of P

    # function DelayLtiSystem(P::StateSpace{Continuous,T, MT}, Tau::Vector{T})
    #     if ControlSystems.noutputs(P) < length(Tau) ||
    #         ControlSystems.noutputs(P) < length(Tau)
    #         error("Length of time-vector is too long given the size of the partitioned system P.")
    #     end
    #     new{T}(P, Tau)
    # end
end

timeevol(sys::DelayLtiSystem) = timeevol(sys.P)

# QUESTION: would psys be a good standard variable name for a PartionedStateSpace
#           and perhaps dsys for a delayed system, (ambigous with discrete system though)
"""
    DelayLtiSystem{T, S}(sys::StateSpace, Tau::AbstractVector{S}=Float64[]) where {T <: Number, S <: Real}

Create a delayed system by speciying both the system and time-delay vector. NOTE: if you want to create a system with simple input or output delays, use the Function `delay(τ)`.
"""
function DelayLtiSystem{T,S}(sys::StateSpace, Tau::AbstractVector{S} = Float64[]) where {T<:Number,S<:Real}
    nu = ninputs(sys) - length(Tau)
    ny = noutputs(sys) - length(Tau)

    if nu < 0  || ny < 0
        throw(ArgumentError("The delay vector of length $length(Tau) is too long."))
    end

    psys = PartionedStateSpace{StateSpace{Continuous,T}}(sys, nu, ny)
    DelayLtiSystem{T,S}(psys, Tau)
end
# For converting DelayLtiSystem{T,S} to different T
DelayLtiSystem{T}(sys::DelayLtiSystem) where {T} = DelayLtiSystem{T}(PartionedStateSpace{StateSpace{Continuous,T}}(sys.P), Float64[])
DelayLtiSystem{T}(sys::StateSpace) where {T} = DelayLtiSystem{T, Float64}(sys, Float64[])

# From StateSpace, infer type
DelayLtiSystem(sys::StateSpace{Continuous,T}, Tau::Vector{S}) where {T,S} = DelayLtiSystem{T,S}(sys, Tau)
DelayLtiSystem(sys::StateSpace{Continuous,T}) where {T} = DelayLtiSystem{T,T}(sys, T[])

# From TransferFunction, infer type TODO Use proper constructor instead of convert here when defined
DelayLtiSystem(sys::TransferFunction{TE,S}) where {TE,T,S<:SisoTf{T}} = DelayLtiSystem{T}(convert(StateSpace{Continuous,T}, sys))


# TODO: Think through these promotions and conversions
Base.promote_rule(::Type{AbstractMatrix{T1}}, ::Type{DelayLtiSystem{T2,S}}) where {T1<:Number,T2<:Number,S} = DelayLtiSystem{promote_type(T1,T2),S}
Base.promote_rule(::Type{T1}, ::Type{DelayLtiSystem{T2,S}}) where {T1<:Number,T2<:Number,S} = DelayLtiSystem{promote_type(T1,T2),S}

Base.promote_rule(::Type{<:StateSpace{<:TimeEvolution,T1}}, ::Type{DelayLtiSystem{T2,S}}) where {T1,T2,S} = DelayLtiSystem{promote_type(T1,T2),S}
Base.promote_rule(::Type{<:TransferFunction}, ::Type{DelayLtiSystem{T,S}}) where {T,S} = DelayLtiSystem{T,S}
#Base.promote_rule(::Type{<:UniformScaling}, ::Type{S}) where {S<:DelayLtiSystem} = DelayLtiSystem{T,S}

function Base.convert(::Type{DelayLtiSystem{T,S}}, sys::StateSpace) where {T,S}
    DelayLtiSystem{T,S}(sys)
end
function Base.convert(::Type{DelayLtiSystem{T1,S}}, d::T2) where {T1,T2 <: Number,S}
    DelayLtiSystem{T1,S}(StateSpace(T1(d)))
end
function Base.convert(::Type{DelayLtiSystem{T1,S}}, d::T2) where {T1,T2 <: AbstractArray,S}
    DelayLtiSystem{T1,S}(StateSpace(T1.(d)))
end

Base.convert(::Type{<:DelayLtiSystem}, sys::TransferFunction)  = DelayLtiSystem(sys)
# Catch convertsion between T
Base.convert(::Type{V}, sys::DelayLtiSystem)  where {T, V<:DelayLtiSystem{T}} =
    sys isa V ? sys : V(StateSpace{Continuous,T}(sys.P.P), sys.Tau)



function *(sys::DelayLtiSystem, n::Number)
    new_C = [sys.P.C1*n; sys.P.C2]
    new_D = [sys.P.D11*n sys.P.D12*n; sys.P.D21 sys.P.D22]
    return DelayLtiSystem(StateSpace(sys.P.A, sys.P.B, new_C, new_D, sys.P.timeevol), sys.Tau)
end
*(n::Number, sys::DelayLtiSystem) = *(sys, n)

function +(sys::DelayLtiSystem{T1,S}, n::T2) where {T1,T2<:Number,S}
    T = promote_type(T1,T2)
    if T == T1 # T2 can be stored in sys
        +(sys, T1(n))
    else # We need to upgrade sys
        +(DelayLtiSystem{T,S}(sys), T(n))
    end
end

# Efficient addition
function +(sys::DelayLtiSystem{T}, n::T) where {T<:Number}
    ny, nu = size(sys)
    ssold = sys.P.P
    # Add to direct term from input to output
    new_D = copy(ssold.D)
    new_D[1:ny, 1:nu] .+= n

    pnew = PartionedStateSpace(StateSpace(ssold.A, ssold.B, ssold.C, new_D, 0.0), ny, nu)
    DelayLtiSystem(pnew, sys.Tau)
end

# Efficient subtraction with number
-(sys::DelayLtiSystem, n::T) where {T <:Number} = +(sys, -n)
-(n::T, sys::DelayLtiSystem) where {T <:Number} = +(-sys, n)
function /(anything, sys::DelayLtiSystem)
    all(iszero, sys.Tau) || error("A delayed system can not be inverted. Consider use of the function `feedback`.")
    /(anything, sys.P.P) # If all delays are zero, invert the inner system
end


# Test equality (of realizations)
function ==(sys1::DelayLtiSystem, sys2::DelayLtiSystem)
    all(getfield(sys1, f) == getfield(sys2, f) for f in fieldnames(DelayLtiSystem))
end


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
        sys.P.timeevol), sys.Tau)
end

function Base.show(io::IO, sys::DelayLtiSystem)
    println(io, typeof(sys))

    print(io, "\nP: ")
    show(io, sys.P.P)

    println(io, "\n\nDelays: $(sys.Tau)")
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

"""
    delay(T::Type{<:Number}, tau)

Create a pure time delay of length `τ` of type `T`.

The type `T` defaults to `promote_type(Float64, typeof(tau))`
"""

function delay(T::Type{<:Number}, τ)
    return DelayLtiSystem(ControlSystems.ss([zero(T) one(T); one(T) zero(T)], Continuous()), [T(τ)])
end
delay(τ::Number) = delay(promote_type(Float64,eltype(τ)), τ)


"""
    exp(G::TransferFunction{Continuous,SisoRational})

Create a time delay of length `tau` with `exp(-τ*s)` where `s=tf("s")` and `τ` > 0.

See also: [`delay`](@ref) which is arguably more conenient than this function.
"""
function Base.exp(G::TransferFunction{Continuous,<:SisoRational})
    if size(G.matrix) != (1,1) && iscontinuous(G)
        error("G must be a continuous-time scalar transfer function. Consider using `delay` instead.")
    end
    G_siso = G.matrix[1,1]

    if !(Polynomials.degree(G_siso.den) == 0 && Polynomials.degree(G_siso.num) == 1
        && G_siso.num[1]/G_siso.den[0] < 0 && G_siso.num[0] == 0)
        error("Input must be of the form -τ*s, τ>0. Consider using `delay` instead.")
    end

    return delay(-G_siso.num[1] / G_siso.den[0])
end
