abstract type LTISystem{T <:TimeType} <: AbstractSystem end
+(sys1::LTISystem, sys2::LTISystem) = +(promote(sys1, sys2)...)
-(sys1::LTISystem, sys2::LTISystem) = -(promote(sys1, sys2)...)
*(sys1::LTISystem, sys2::LTISystem) = *(promote(sys1, sys2)...)
/(sys1::LTISystem, sys2::LTISystem) = /(promote(sys1, sys2)...)

# Fallback number
+(sys1::LTISystem, sys2::Number) = +(promote(sys1, sys2)...)
-(sys1::LTISystem, sys2::Number) = -(promote(sys1, sys2)...)
*(sys1::LTISystem, sys2::Number) = *(promote(sys1, sys2)...)
/(sys1::LTISystem, sys2::Number) = /(promote(sys1, sys2)...)

+(sys1::Number, sys2::LTISystem) = +(promote(sys1, sys2)...)
-(sys1::Number, sys2::LTISystem) = -(promote(sys1, sys2)...)
*(sys1::Number, sys2::LTISystem) = *(promote(sys1, sys2)...)
/(sys1::Number, sys2::LTISystem) = /(promote(sys1, sys2)...)

# Fallback Matrix
+(sys1::LTISystem, sys2::AbstractMatrix) = +(promote(sys1, sys2)...)
-(sys1::LTISystem, sys2::AbstractMatrix) = -(promote(sys1, sys2)...)
*(sys1::LTISystem, sys2::AbstractMatrix) = *(promote(sys1, sys2)...)
/(sys1::LTISystem, sys2::AbstractMatrix) = /(promote(sys1, sys2)...)

+(sys1::AbstractMatrix, sys2::LTISystem) = +(promote(sys1, sys2)...)
-(sys1::AbstractMatrix, sys2::LTISystem) = -(promote(sys1, sys2)...)
*(sys1::AbstractMatrix, sys2::LTISystem) = *(promote(sys1, sys2)...)
/(sys1::AbstractMatrix, sys2::LTISystem) = /(promote(sys1, sys2)...)

# TODO We should have proper fallbacks for matrices too

feedback(sys1::Union{LTISystem,Number,AbstractMatrix{<:Number}},
         sys2::Union{LTISystem,Number,AbstractMatrix{<:Number}}; kwargs...) = feedback(promote(sys1, sys2)...; kwargs...)

Base.inv(G::LTISystem) = 1/G


"""`issiso(sys)`

Returns `true` if `sys` is SISO, else returns `false`."""
function issiso(sys::LTISystem)
    return ninputs(sys) == 1 && noutputs(sys) == 1
end


"""`iscontinuous(sys)`

Returns `true` if `sys` is continuous, else returns `false`."""
iscontinuous(::LTISystem) = false
iscontinuous(::LTISystem{<:Continuous}) = true
"""`isdiscrete(sys)`

Returns `true` if `sys` is discrete, else returns `false`."""
isdiscrete(::LTISystem) = false
isdiscrete(::LTISystem{<:Discrete}) = true
"""`isstatic(sys)`

Returns `true` if `sys` is static, else returns `false`."""

"""`sampletime(sys)`

Returns the sampletime of a discrete time system, throws error if the system is continuous time or static.
Always ensure `isdiscrete` before using."""
sampletime(sys::LTISystem{<:Discrete}) = sys.time.Ts
sampletime(::LTISystem{<:Continuous}) = error("Continuous system has no sample time")

"""`sampletype(sys)`
Get the sampletype of system. Usually typeof(sys.time)."""
sampletype(sys) = typeof(sys.time) # QUESTION: Perhaps timetype or sampletimetype instead?

common_sample_time(systems::LTISystem...) = common_sample_time(s.time for s in systems)


"""`isstable(sys)`

Returns `true` if `sys` is stable, else returns `false`."""
function isstable(sys::LTISystem)
    if iscontinuous(sys)
        if any(real.(pole(sys)).>=0)
            return false
        end
    else
        if any(abs.(pole(sys)).>=1)
            return false
        end
    end
    return true
end

# Fallback since LTISystem not AbstractArray
Base.size(sys::LTISystem, i::Integer) = size(sys)[i]
