const nonlinear_warning = """
NOTE: The nonlinear functionality in ControlSystems.jl is currently experimental and subject to breaking changes not respecting semantic versioning. Use at your own risk.
"""

"""
    struct HammersteinWienerSystem{T} <: LTISystem

Represents an LTISystem with element-wise static nonlinearities at the inputs and outputs. See `?nonlinearity` for a convenience constructor.
$nonlinear_warning

```
      ┌─────────┐
 y◄───┤         │◄────u
      │    P    │
Δy┌───┤         │◄───┐Δu
  │   └─────────┘    │
  │      ┌───┐       │
  └─────►│ f ├───────┘
         └───┘
```
"""
struct HammersteinWienerSystem{T} <: LFTSystem{Continuous, T}
    P::PartitionedStateSpace{Continuous, StateSpace{Continuous,T}}
    f::Vector{Function} # The length of the vector f implicitly defines the partitionging of P
end

feedback_channel(sys::HammersteinWienerSystem) = sys.f

timeevol(sys::HammersteinWienerSystem) = timeevol(sys.P)

"""
    HammersteinWienerSystem{T, S}(sys::StateSpace, f::Vector{Function}=[]) where {T <: Number}

Create a nonlinear system by speciying both the system and nonlinearity. Users should prefer to use the function [`nonlinearity`](@ref).
"""
function HammersteinWienerSystem{T}(sys::StateSpace, f::Vector{Function} = Function[]) where {T<:Number}
    nu = ninputs(sys) - length(f)
    ny = noutputs(sys) - length(f)

    (nu < 0  || ny < 0) && throw(ArgumentError("The nonlinearity vector of length $(length(f)) is too long."))

    psys = PartitionedStateSpace{Continuous, StateSpace{Continuous,T}}(sys, nu, ny)
    HammersteinWienerSystem{T}(psys, f)
end
# For converting HammersteinWienerSystem{T,S} to different T
HammersteinWienerSystem{T}(sys::HammersteinWienerSystem) where {T} = HammersteinWienerSystem{T}(PartitionedStateSpace{Continuous, StateSpace{Continuous,T}}(sys.P), Function[])
HammersteinWienerSystem{T}(sys::StateSpace) where {T} = HammersteinWienerSystem{T}(sys, Function[])
HammersteinWienerSystem(sys::HammersteinWienerSystem) = sys
# From StateSpace, infer type
HammersteinWienerSystem(sys::StateSpace{Continuous,T}, f::Vector{Function} = Function[]) where {T} = HammersteinWienerSystem{T}(sys, f)

# From TransferFunction, infer type TODO Use proper constructor instead of convert here when defined
HammersteinWienerSystem(sys::TransferFunction{TE,S}) where {TE,T,S<:SisoTf{T}} = HammersteinWienerSystem{T}(convert(StateSpace{Continuous,T}, sys))


Base.promote_rule(::Type{AbstractMatrix{T1}}, ::Type{HammersteinWienerSystem{T2}}) where {T1<:Number,T2<:Number} = HammersteinWienerSystem{promote_type(T1,T2)}
Base.promote_rule(::Type{T1}, ::Type{HammersteinWienerSystem{T2}}) where {T1<:Number,T2<:Number} = HammersteinWienerSystem{promote_type(T1,T2)}

Base.promote_rule(::Type{<:StateSpace{<:TimeEvolution,T1}}, ::Type{HammersteinWienerSystem{T2}}) where {T1,T2} = HammersteinWienerSystem{promote_type(T1,T2)}
Base.promote_rule(::Type{<:TransferFunction{<:Any, ST}}, ::Type{HammersteinWienerSystem{T}}) where {T, ST} = HammersteinWienerSystem{promote_type(T, numeric_type(ST))}

function Base.convert(::Type{HammersteinWienerSystem{T}}, sys::StateSpace) where {T}
    HammersteinWienerSystem{T}(sys)
end
function Base.convert(::Type{HammersteinWienerSystem{T1}}, d::T2) where {T1,T2 <: Number}
    HammersteinWienerSystem{T1}(StateSpace(T1(d)))
end
function Base.convert(::Type{HammersteinWienerSystem{T1}}, d::T2) where {T1,T2 <: AbstractArray}
    HammersteinWienerSystem{T1}(StateSpace(T1.(d)))
end

function Base.convert(::Type{HammersteinWienerSystem{T}}, sys::TransferFunction{TE}) where {T,TE}
    HammersteinWienerSystem{T}(convert(StateSpace{TE, T}, sys))
end
# Catch convertsion between T
Base.convert(::Type{V}, sys::HammersteinWienerSystem)  where {T, V<:HammersteinWienerSystem{T}} =
    sys isa V ? sys : V(StateSpace{Continuous,T}(sys.P.P), sys.f)

function Base.show(io::IO, sys::HammersteinWienerSystem)
    println(io, typeof(sys))

    print(io, "\nP: ")
    show(io, sys.P.P)

    print(io, "\n\nNonlinearities: $(sys.f)")
end


function Base.getproperty(sys::HammersteinWienerSystem, s::Symbol)
    s ∈ fieldnames(typeof(sys)) && return getfield(sys, s)
    s ∈ propertynames(getfield(sys, :P)) && return getproperty(getfield(sys, :P), s)
    throw(ArgumentError("$(typeof(sys)) has no property named $s"))
end


"""
    nonlinearity(T, f)

Create a pure nonlinearity.

$nonlinear_warning

The type `T` defaults to `promote_type(Float64, typeof(f))`

# Example:
Create a LTI system with an input nonlinearity of `f`
```julia
tf(1, [1, 1])*nonlinearity(x->clamp(x, -1, 1))
```

See also predefined nonlinearities [`saturation`](@ref), [`offser`](@ref).
"""
function nonlinearity(::Type{T}, f) where T <: Number
    return HammersteinWienerSystem(ControlSystems.ss([zero(T) one(T); one(T) zero(T)], Continuous()), Function[f])
end

function nonlinearity(::Type{T}, fs::AbstractVector) where T <: Number
    return append(nonlinearity.(T, fs)...)
end

nonlinearity(f) = nonlinearity(Float64, f)