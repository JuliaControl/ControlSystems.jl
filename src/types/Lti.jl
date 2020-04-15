abstract type LTISystem <: AbstractSystem end
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


"""`is_continuous_time(sys)`

Returns `true` for a continuous-time system `sys`, else returns `false`."""
is_continuous_time(sys::LTISystem) = sampletime(sys) isa Continuous
"""`is_discrete_time(sys)`

Returns `true` for a discrete-time system `sys`, else returns `false`."""
is_discrete_time(sys::LTISystem) = sampletime(sys) isa Discrete


function Base.getproperty(sys::LTISystem, s::Symbol)
    if s === :Ts
        # if !is_discrete_time(sys) # NOTE this line seems to be breaking inference of is_discrete_time (is there a test for this?)
        if !is_discrete_time(sys)
            @warn "Getting sampletime 0.0 for non-discrete systems is deprecated. Check `is_discrete_time` before trying to access sampletime."
            return 0.0
        else
            return sampletime(sys).Ts
        end
    else
        return getfield(sys, s)
    end
end


"""`timetype(sys)`
Get the timetype of system. Usually typeof(sys.sampletime)."""
timetype(sys) = typeof(sys.sampletime)

common_sampletime(systems::LTISystem...) = common_sampletime(sampletime(sys) for sys in systems)


"""`isstable(sys)`

Returns `true` if `sys` is stable, else returns `false`."""
function isstable(sys::LTISystem)
    if is_continuous_time(sys)
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
