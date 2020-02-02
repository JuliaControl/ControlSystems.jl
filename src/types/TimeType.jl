
abstract type TimeType end

const UNDEF_SAMPLETIME = -1

struct Discrete{T} <: TimeType
    Ts::T
    function Discrete{T}(Ts::T) where T
        if Ts <= 0 && Ts != UNDEF_SAMPLETIME
            throw(ErrorException("Creating a continuous time system by setting sample time to 0 is no longer supported."))
        end
        new{T}(Ts)
    end
end
Discrete{T}(x) where T = Discrete{T}(T(x))

struct Continuous <: TimeType end
struct Static <: TimeType end
# Basic Constructors
Discrete(x::T) where T<:Number = Discrete{T}(x)
Continuous(x::Continuous) = x
Static(x::Static) = x
# Simple converseion
Discrete{T}(x::Discrete) where T = Discrete{T}(x.Ts)
Discrete{T}(::Static) where T = Discrete{T}(UNDEF_SAMPLETIME)
Continuous(::Static) = Continuous()
# Default to checing type
iscontinuous(x::TimeType) = x isa Continuous
isdiscrete(x::TimeType) = x isa Discrete
isstatic(x::TimeType) = x isa Static

# Interface for getting sample time
sampletime(x::Discrete) = x.Ts
sampletime(x::Continuous) = error("Continuous system has no sample time")
sampletime(x::Static) = error("Static system has no sample time")

unknown_time(::Type{Discrete{T}}) where T = Discrete{T}(UNDEF_SAMPLETIME)
unknown_time(::Type{Continuous}) where T = Continuous()
unknown_time(::Type{Static}) where T = Static()

# Promotion
# Fallback to give useful error
Base.promote_rule(::Type{<:Discrete}, ::Type{Continuous}) = throw(ErrorException("Sampling time mismatch"))
Base.promote_rule(::Type{T}, ::Type{Static}) where T <: TimeType = T
Base.promote_rule(::Type{Discrete{T1}}, ::Type{Discrete{T2}}) where {T1,T2}= Discrete{promote_type(T1,T2)}

Base.convert(::Type{Discrete{T1}}, Ts::Discrete{T2}) where {T1,T2} = Discrete{T1}(Ts.Ts)

# Promoting two or more systems systems should promote sample times
common_sample_time(x::TimeType) = x
common_sample_time(x::TimeType, y::TimeType) = throw(ErrorException("Sampling time mismatch"))
common_sample_time(x::TimeType, y::TimeType, z...) = common_sample_time(common_sample_time(x, y), z...)
common_sample_time(a::Base.Generator) = reduce(common_sample_time, a)

function common_sample_time(x::Discrete{T1}, y::Discrete{T2}) where {T1,T2}
    if x != y && x.Ts != UNDEF_SAMPLETIME && y.Ts != UNDEF_SAMPLETIME
        throw(ErrorException("Sampling time mismatch"))
    end

    if x.Ts == UNDEF_SAMPLETIME
        return Discrete{promote_type(T1,T2)}(y)
    else
        return Discrete{promote_type(T1,T2)}(x)
    end
end

common_sample_time(x::Continuous, ys::Continuous...) = Continuous()
common_sample_time(x::Static, ys::Static...) = Static()

# Combinations with static
common_sample_time(x::TimeType, ys::Static) = x
common_sample_time(x::Static, y::TimeType) = y

# Check equality and similarity
==(x::TimeType, y::TimeType) = false
==(x::Discrete, y::Discrete) = (x.Ts == y.Ts)
==(::Continuous, ::Continuous) = true
==(::Static, ::Static) = true

# Equality with static is true
==(::Static, ::TimeType) = true
==(::TimeType, ::Static) = true

isapprox(x::TimeType, y::TimeType, args...; kwargs...) = false
isapprox(x::Discrete, y::Discrete, args...; kwargs...) = isapprox(x.Ts, y.Ts, args...; kwargs...)
isapprox(::Continuous, ::Continuous, args...; kwargs...) = true
isapprox(::Static, ::Static, args...; kwargs...) = true

# Equality with static is true
isapprox(::Static, ::TimeType, args...; kwargs...) = true
isapprox(::TimeType, ::Static, args...; kwargs...) = true
