module ReferenceGenerators

using Distributions

export PWCGenerator, 
        ChirpGenerator, 
        RSGenerator, 
        SawtoothGenerator, 
        SquareGenerator,
        GNGenerator,
        NoiseGenerator,
        RGCollector

#PWC
"""
r = PWCGenerator(t::T, v::T)

Generates a PIECE-WISE CONSTANT reference generator where a call
    ref = r(t)
returns the reference value at time t.

NOTE: Input vectors has to be of equal length
INPUT:
    t::AbstractVector{<:Real} = the time points in which the constant value changes
    v::AbstractVector{<:Real} = the constant value the reference changes to
"""
struct PWCGenerator{T <: AbstractVector{<:Real}} <: Function
    t::T
    v::T
end

function (ref::PWCGenerator)(t::Real)
    if t < ref.t[1]
        return 0
    elseif t >= ref.t[end]
        return ref.v[end]
    end

    for (idx, time) in enumerate(ref.t)
        if time > t && idx != 1
            return ref.v[idx - 1]
        end
    end
end

#Chirp
"""
r = ChirpGenerator(f::T = t -> t)

Generates a CHIRP signal generator where a call
    ref = r(t)
returns the chirp value at time t.

DEFAULT: 
    Linear frequency chirp, f(t) = t

INPUT:
    f::Function = the frequency function of the chirp signal (f = function of 1 variable)
"""
struct ChirpGenerator{T <: Function} <: Function
    f::T
end

ChirpGenerator() = ChirpGenerator(t -> t)

function (ref::ChirpGenerator)(t::Real)
    return sin(ref.f(t)*t)
end

#RS
"""
r = RSGenerator(fun::T, f::S)

Generates a REPEATING SEQUENCE generator where a call
    ref = r(t)
returns the repeating sequence's value at time t.

INPUT:
    fun::Function = the function to repeat 
    f::Real       = The frequency of which to repeat the function. If f = 2 the
                    RSGenerator will repeat the first 0.5 time units of fun infinitely
"""
struct RSGenerator{T <: Function, S <: Real} <: Function
    fun::T
    f::S
end

function (ref::RSGenerator)(t::Real)
    # Reduces the time till it lies within the interval 0 <= t <= ref.t
    T = 1 / ref.f
    t = mod(t, T)

    return ref.fun(t)
end

#Sawtooth
"""
r = SawtoothGenerator(f::T, p::T = 0, frac::Function = t -> t - floor(t))

Generates a SAWTOOTH signal generator where a call
    ref = r(t)
returns the sawtooth value at time t.

DEFAULT:
    p::T = The phase of the sawtooth signal is default set to 0
    frac::Function = The sawtooth function. Please do not alter this

INPUT:
    f::T = The frequency of which to repeat the sawtooth signal
    p::T = The phase of the sawtooth signal
    frac::Function = The sawtooth function. Please do not alter this
"""
struct SawtoothGenerator{T <: Real} <: Function
    f::T #Frequency
    p::T #Phase
end

SawtoothGenerator(f)    = SawtoothGenerator(f, 0)
SawtoothGenerator(f, p) = SawtoothGenerator(f, p)

function (ref::SawtoothGenerator)(t::Real)
    return (ref.f*t+ref.p) - floor(ref.f*t+ref.p)
end


#Square
"""
r = SquareGenerator(f::T, A::T = 1)

Generates a SQUARE signal generator where a call
    ref = r(t)
returns the square value at time t.

DEFAULT:
    A::T = The amplitude of the square wave is initially set to 1

INPUT:
    f::T = The frequency of the square wave
    A::T = The amplitude of the square wave
"""
struct SquareGenerator{T <: Real} <: Function
    f::T #Frequency
    A::T #Amplitude
end

SquareGenerator(f) = SquareGenerator(f, (1)::typeof(f))

function (ref::SquareGenerator)(t::Real)
    T = 1 / ref.f
    t = mod(t, T)

    return (t < T/2) ? ref.A : -ref.A
end

#GN
"""
r = GNGenerator(v::T, m::T = 0, N = Normal(m, v))

Generates a GAUSSIAN NOISE signal generator where a call
    ref = r(t)
returns the gaussian noise value at time t.

DEFAULT:
    m::T = The mean is default set to 0
    N = N is initially (and should always) be set to to the normal function
        from the Distributions package 

INPUT:
    v::T = The variance of the GNGenerator
    m::T = The mean of the GNGenerator
    N    = The normal distribution
"""
struct GNGenerator{T <: Real} <: Function
    v::T #Variance
    m::T #Mean
    N::Distribution    #TODO: How to fix this?
end

GNGenerator(v, m) = GNGenerator(m, v, Normal(m, v))
GNGenerator(v)    = GNGenerator((0)::typeof(v), v, Normal((0)::typeof(v), v))

function (ref::GNGenerator)(t::Real)
    return rand(ref.N)
end

#Noise
"""
`r = NoiseGenerator(rng::T)`

Generates a NOISE signal generator where a call
    ref = r(t)
returns the given noise value at time t.

INPUT:
    rng::Distribution = The distribution to sample the noise from
"""
struct NoiseGenerator{T <: Distribution} <: Function
    rng::T
end

function (ref::NoiseGenerator)(t::Real)
    return rand(ref.rng)
end

#Collector
"""
r = RGCollector(generators::T)

Generates a REFERENCE GENERATOR COLLECTOR where a call
    ref = r(t)
returns the reference generators' value at time t.

INPUT:
    generators::AbstractVector{Function} = A vector of reference generators
"""
struct RGCollector{T <: Tuple{Vararg{Function}}} <: Function
    generators::T
end

RGCollector(args::Function...) = RGCollector((args...,))

function (rgc::RGCollector)(t::Real)
    return [generator(t) for generator in rgc.generators]
end

end #module
