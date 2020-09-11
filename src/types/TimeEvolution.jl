
abstract type TimeEvolution end

const UNDEF_SAMPLEPETIME = -1 # For handling promotion of Matrix to LTISystem

struct Discrete{T} <: TimeEvolution
    Ts::T
    function Discrete{T}(Ts::T) where T
        if Ts <= 0 && Ts != UNDEF_SAMPLEPETIME
            throw(ErrorException("Creating a continuous time system by setting sample time to 0 is no longer supported."))
        end
        new{T}(Ts)
    end
end
Discrete{T}(x) where T = Discrete{T}(T(x))

struct Continuous <: TimeEvolution end
# Basic Constructors
Discrete(x::T) where T<:Number = Discrete{T}(x)
Continuous(x::Continuous) = x
# Simple converseion
Discrete{T}(x::Discrete) where T = Discrete{T}(x.Ts)


undef_sampletime(::Type{Discrete{T}}) where T = Discrete{T}(UNDEF_SAMPLEPETIME)
undef_sampletime(::Type{Continuous}) where T = Continuous()


# Promotion
# Fallback to give useful error
Base.promote_rule(::Type{<:Discrete}, ::Type{Continuous}) = throw(ErrorException("Sampling time mismatch"))
Base.promote_rule(::Type{Discrete{T1}}, ::Type{Discrete{T2}}) where {T1,T2}= Discrete{promote_type(T1,T2)}

Base.convert(::Type{Discrete{T1}}, x::Discrete{T2}) where {T1,T2} = Discrete{T1}(x.Ts)

# Promoting two or more systems systems should promote sample times
common_timeevol(x::TimeEvolution) = x
common_timeevol(x::TimeEvolution, y::TimeEvolution) = throw(ErrorException("Sampling time mismatch"))
common_timeevol(x::TimeEvolution, y::TimeEvolution, z...) = common_timeevol(common_timeevol(x, y), z...)
common_timeevol(a::Base.Generator) = reduce(common_timeevol, a)

function common_timeevol(x::Discrete{T1}, y::Discrete{T2}) where {T1,T2}
    if x != y && x.Ts != UNDEF_SAMPLEPETIME && y.Ts != UNDEF_SAMPLEPETIME
         throw(ErrorException("Sampling time mismatch"))
    end

    if x.Ts == UNDEF_SAMPLEPETIME
        return Discrete{promote_type(T1,T2)}(y)
    else
        return Discrete{promote_type(T1,T2)}(x)
    end
end

common_timeevol(x::Continuous, ys::Continuous...) = Continuous()

# Check equality
==(x::TimeEvolution, y::TimeEvolution) = false
==(x::Discrete, y::Discrete) = (x.Ts == y.Ts)
==(::Continuous, ::Continuous) = true

isapprox(x::TimeEvolution, y::TimeEvolution, args...; kwargs...) = false
isapprox(x::Discrete, y::Discrete, args...; kwargs...) = isapprox(x.Ts, y.Ts, args...; kwargs...)
isapprox(::Continuous, ::Continuous, args...; kwargs...) = true
