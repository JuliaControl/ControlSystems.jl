
"""
    struct DelayLtiSystem{T, S <: Real} <: LTISystem

Represents an LTISystem with internal time-delay. See `?delay` for a convenience constructor.
"""
struct DelayLtiSystem{T,S<:Real} <: LFTSystem{Continuous, T}
    P::PartitionedStateSpace{Continuous, StateSpace{Continuous,T}}
    Tau::Vector{S} # The length of the vector tau implicitly defines the partitionging of P
end

feedback_channel(sys::DelayLtiSystem) = sys.Tau


"""
    DelayLtiSystem{T, S}(sys::StateSpace, Tau::AbstractVector{S}=Float64[]) where {T <: Number, S <: Real}

Create a delayed system by specifying both the system and time-delay vector. NOTE: if you want to create a system with simple input or output delays, use the Function `delay(τ)`.
"""
function DelayLtiSystem{T,S}(sys::StateSpace, Tau::AbstractVector{S} = Float64[]) where {T<:Number,S<:Real}
    nu = ninputs(sys) - length(Tau)
    ny = noutputs(sys) - length(Tau)

    if nu < 0  || ny < 0
        throw(ArgumentError("The delay vector of length $length(Tau) is too long."))
    end
    csys = convert(StateSpace{Continuous,T}, sys)
    psys = PartitionedStateSpace{Continuous, StateSpace{Continuous,T}}(csys, nu, ny)
    DelayLtiSystem{T,S}(psys, Tau)
end
# For converting DelayLtiSystem{T,S} to different T
DelayLtiSystem{T}(sys::DelayLtiSystem) where {T} = DelayLtiSystem{T}(PartitionedStateSpace{Continuous, StateSpace{Continuous,T}}(sys.P), Float64[])
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
Base.promote_rule(::Type{<:TransferFunction{<:Any, ST}}, ::Type{DelayLtiSystem{T,S}}) where {T,S, ST} = DelayLtiSystem{promote_type(T, numeric_type(ST)),S}
#Base.promote_rule(::Type{<:UniformScaling}, ::Type{S}) where {S<:DelayLtiSystem} = DelayLtiSystem{T,S}

function Base.convert(::Type{DelayLtiSystem{T,S}}, sys::StateSpace) where {T,S}
    DelayLtiSystem{T,S}(convert(StateSpace{Continuous,T}, sys))
end
function Base.convert(::Type{DelayLtiSystem{T1,S}}, d::T2) where {T1,T2 <: Number,S}
    DelayLtiSystem{T1,S}(StateSpace(T1(d)))
end
function Base.convert(::Type{DelayLtiSystem{T1,S}}, d::T2) where {T1,T2 <: AbstractArray,S}
    DelayLtiSystem{T1,S}(StateSpace(T1.(d)))
end

function Base.convert(::Type{DelayLtiSystem{T,S}}, sys::TransferFunction{TE}) where {T,S,TE}
    if issiso(sys) && length(numvec(sys.matrix[1,1])) > length(denvec(sys.matrix[1,1]))
        error("The transfer function is not proper and can not be converted to a DelayLtiSystem type. If you tried to form the system `exp(sL) * B / A` where `B / A` is proper, add parenthesis to make it `exp(sL) * (B / A)`.")
    end
       
    DelayLtiSystem{T,S}(convert(StateSpace{TE, T}, sys))
end
# Catch conversion between T
Base.convert(::Type{V}, sys::DelayLtiSystem)  where {T, V<:DelayLtiSystem{T}} =
    sys isa V ? sys : V(StateSpace{Continuous,T}(sys.P.P), sys.Tau)



# Efficient addition
function +(sys::DelayLtiSystem{<:Any,T}, n::T) where {T<:Number}
    ny, nu = size(sys)
    ssold = sys.P.P
    # Add to direct term from input to output
    new_D = copy(ssold.D)
    new_D[1:ny, 1:nu] .+= n

    pnew = PartitionedStateSpace(StateSpace(ssold.A, ssold.B, ssold.C, new_D, Continuous()), ny, nu)
    DelayLtiSystem(pnew, feedback_channel(sys))
end


function +(sys::DelayLtiSystem{T1,S}, n::T2) where {T1,T2<:Number,S}
    T = promote_type(T1,T2)
    if T == T1 # T2 can be stored in sys
        +(sys, T1(n))
    else # We need to upgrade sys
        +(DelayLtiSystem{T,S}(sys), T(n))
    end
end


"""
    /(G1, G2::DelayLtiSystem)

Compute ``G_1 * G_2^{-1} where ``G_2`` is a DelayLtiSystem.
Throws a SingularException if ``G_2`` is not invertible.
"""
function /(anything, sys::DelayLtiSystem)
    ny,nu = size(sys)
    ny == nu || error("The denominator system must be square")
    return anything * feedback(I(nu), sys - I(nu))
end

for other_type in [:Number, :AbstractMatrix, :LTISystem]
    @eval /(a::$other_type, sys::DelayLtiSystem) = invoke(/, Tuple{Any, DelayLtiSystem}, a, sys)
end


function Base.show(io::IO, sys::DelayLtiSystem)
    println(io, typeof(sys))

    print(io, "\nP: ")
    show(io, sys.P.P)

    print(io, "\n\nDelays: $(sys.Tau)")
end



"""
    delay(tau)
    delay(tau, Ts)
    delay(T::Type{<:Number}, tau)
    delay(T::Type{<:Number}, tau, Ts)

Create a pure time delay of length `τ` of type `T`.

The type `T` defaults to `promote_type(Float64, typeof(tau))`.

If `Ts` is given, the delay is discretized with sampling time `Ts` and a discrete-time StateSpace object is returned.

# Example:
Create a LTI system with an input delay of `L`
```julia
L = 1
tf(1, [1, 1])*delay(L)
s = tf("s")
tf(1, [1, 1])*exp(-s*L) # Equivalent to the version above
```
"""
function delay(T::Type{<:Number}, τ)
    return DelayLtiSystem(ControlSystemsBase.ss([zero(T) one(T); one(T) zero(T)], Continuous()), [T(τ)])
end
delay(τ::Number) = delay(promote_type(Float64,eltype(τ)), τ)

delay(τ::Number, Ts::Number) = c2d(delay(τ), Ts)
delay(T::Type{<:Number}, τ::Number, Ts::Number) = c2d(delay(T, τ), Ts)


"""
    exp(G::TransferFunction{Continuous,SisoRational})

Create a time delay of length `tau` with `exp(-τ*s)` where `s=tf("s")` and `τ` > 0.

See also: [`delay`](@ref) which is arguably more convenient than this function.
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
