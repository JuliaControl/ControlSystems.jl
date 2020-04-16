
abstract type TimeType end

const UNDEF_SAMPLETIME = -1 # For handling promotion of Matrix to LTISystem

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
# Basic Constructors
Discrete(x::T) where T<:Number = Discrete{T}(x)
Continuous(x::Continuous) = x
# Simple converseion
Discrete{T}(x::Discrete) where T = Discrete{T}(x.Ts)

undef_sampletime(::Type{Discrete{T}}) where T = Discrete{T}(UNDEF_SAMPLETIME)
undef_sampletime(::Type{Continuous}) where T = Continuous()

sampletime(x::Discrete) = x.Ts
sampletime(x::Continuous) = error("Continuous system has no sample time")

# Promotion
# Fallback to give useful error
Base.promote_rule(::Type{<:Discrete}, ::Type{Continuous}) = throw(ErrorException("Sampling time mismatch"))
Base.promote_rule(::Type{Discrete{T1}}, ::Type{Discrete{T2}}) where {T1,T2}= Discrete{promote_type(T1,T2)}

Base.convert(::Type{Discrete{T1}}, x::Discrete{T2}) where {T1,T2} = Discrete{T1}(x.Ts)

# Promoting two or more systems systems should promote sample times
common_time(x::TimeType) = x
common_time(x::TimeType, y::TimeType) = throw(ErrorException("Sampling time mismatch"))
common_time(x::TimeType, y::TimeType, z...) = common_time(common_time(x, y), z...)
common_time(a::Base.Generator) = reduce(common_time, a)

function common_time(x::Discrete{T1}, y::Discrete{T2}) where {T1,T2}
    if x != y && x.Ts != UNDEF_SAMPLETIME && y.Ts != UNDEF_SAMPLETIME
         throw(ErrorException("Sampling time mismatch"))
    end

    if x.Ts == UNDEF_SAMPLETIME
        return Discrete{promote_type(T1,T2)}(y)
    else
        return Discrete{promote_type(T1,T2)}(x)
    end
end

common_time(x::Continuous, ys::Continuous...) = Continuous()

# Check equality
==(x::TimeType, y::TimeType) = false
==(x::Discrete, y::Discrete) = (x.Ts == y.Ts)
==(::Continuous, ::Continuous) = true

isapprox(x::TimeType, y::TimeType, args...; kwargs...) = false
isapprox(x::Discrete, y::Discrete, args...; kwargs...) = isapprox(x.Ts, y.Ts, args...; kwargs...)
isapprox(::Continuous, ::Continuous, args...; kwargs...) = true
