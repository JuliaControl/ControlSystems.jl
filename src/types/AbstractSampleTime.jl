
abstract type AbstractSampleTime end

struct Discrete{T} <: AbstractSampleTime
    Ts::T
    function Discrete{T}(Ts::T) where T
        if Ts <= 0
            throw(ArgumentError("Creating a continuous time system by setting sample time to 0 is no longer supported."))
        end
        new{T}(Ts)
    end
end
struct Continuous <: AbstractSampleTime end
struct Static <: AbstractSampleTime end
# Basic Constructors
Discrete(x::T) where T<:Number = Discrete{T}(x)
Continuous(x::Continuous) = x
Static(x::Static) = x
# Simple converseion
Discrete{T}(x::Discrete) where T = Discrete{T}(x.Ts)
# Check time sampling by type
iscontinuous(::Type{<:AbstractSampleTime}) = false
iscontinuous(::Type{Continuous}) = true
isdiscrete(::Type{<:AbstractSampleTime}) = false
isdiscrete(::Type{<:Discrete}) = true
isstatic(::Type{<:AbstractSampleTime}) = false
isstatic(::Type{Static}) = true
# Default to checing type
iscontinuous(x::AbstractSampleTime) = iscontinuous(typeof(x))
isdiscrete(x::AbstractSampleTime) = isdiscrete(typeof(x))
isstatic(x::AbstractSampleTime) = isstatic(typeof(x))

# Interface for getting sample time
sampletime(x::Discrete) = x.Ts
sampletime(x::Continuous) = error("Continuous system has no sample time")
sampletime(x::Static) = error("Static system has no sample time")

# Promotion
Base.promote_rule(::Type{Continuous}, ::Type{Continuous}) = Continuous
Base.promote_rule(::Type{Discrete{T1}}, ::Type{Discrete{T2}}) where {T1,T2}= Discrete{promote_type(T1,T2)}

Base.convert(::Type{Discrete{T1}}, Ts::Discrete{T2}) where {T1,T2} = Discrete{T1}(Ts.Ts)

# Promoting two or more systems systems should promote sample times
ts_same(x::AbstractSampleTime, y::AbstractSampleTime) = throw(ArgumentError("Sampling time mismatch"))
ts_same(x::AbstractSampleTime, y::AbstractSampleTime, z...) = ts_same(ts_same(x, y), z...)

function ts_same(x::Discrete{T}, y::Discrete{S}) where {T,S}
    if x != y
        throw(ArgumentError("Sampling time mismatch"))
    end
    return Discrete{promote_type(T,S)}(x)
end
ts_same(x::Continuous, y::Continuous) = Continuous()



# Check equality and similarity
==(x::AbstractSampleTime, y::AbstractSampleTime) = false
==(x::Discrete, y::Discrete) = (x.Ts == y.Ts)
==(::Continuous, ::Continuous) = true
==(::Static, ::Static) = true

≈(x::AbstractSampleTime, y::AbstractSampleTime) = false
≈(x::Discrete, y::Discrete) = (x.Ts ≈ y.Ts)
≈(::Continuous, ::Continuous) = true
≈(::Static, ::Static) = true
