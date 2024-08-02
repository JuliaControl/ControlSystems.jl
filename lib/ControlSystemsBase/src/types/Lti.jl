abstract type LTISystem{TE<:TimeEvolution} <: AbstractSystem end
+(sys1::LTISystem, sys2::LTISystem) = +(promote(sys1, sys2)...)
-(sys1::LTISystem, sys2::LTISystem) = -(promote(sys1, sys2)...)
function *(sys1::LTISystem, sys2::LTISystem)
    zeroexcess(sys) = length(numvec(sys)[]) - length(denvec(sys)[])
    try
        *(promote(sys1, sys2)...)
    catch e
        e isa ImproperException || rethrow()
        if sys1 isa AbstractStateSpace && sys2 isa TransferFunction
            issiso(sys2) || rethrow() # Can't invert MINO tf
            zeroexcess(sys2) <= sys1.nx || rethrow() # Quotient can't be proper
            Q = sys1 / ss(inv(sys2))
        elseif sys1 isa TransferFunction && sys2 isa AbstractStateSpace
            issiso(sys1) || rethrow() # Can't invert MINO tf
            zeroexcess(sys1) <= sys2.nx || rethrow() # Quotient can't be proper
            Q = ss(inv(sys1)) \ sys2
        else
            rethrow()
        end
        Q = balance_statespace(Q)[1]
        if max(maximum(abs, Q.A), maximum(abs, Q.B), maximum(abs, Q.C), maximum(abs, Q.D)) > 1e7
            @warn "Possible numerical instability detected: Multiplication of a statespace system and a non-proper transfer function may result in numerical inaccuracy. Verify result carefully, and consider making use of DescriptorSystems.jl to represent this product as a DescriptorSystem with non-unit descriptor matrix if result is inaccurate."
        end
        return Q
    end
end
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

Base.:^(sys::LTISystem, p::Integer) = Base.power_by_squaring(sys, p)

# TODO We should have proper fallbacks for matrices too

feedback(sys1::Union{LTISystem,Number,AbstractMatrix{<:Number}},
         sys2::Union{LTISystem,Number,AbstractMatrix{<:Number}}; kwargs...) = feedback(promote(sys1, sys2)...; kwargs...)

Base.inv(G::LTISystem) = 1/G


"""`issiso(sys)`

Returns `true` if `sys` is SISO, else returns `false`."""
function issiso(sys::LTISystem)
    return ninputs(sys) == 1 && noutputs(sys) == 1
end

"""`timeevol(sys::LTISystem)`
Get the timeevol::TimeEvolution from system `sys`, usually sys.timeevol """
timeevol(sys::LTISystem) = sys.timeevol

"""`iscontinuous(sys)`

Returns `true` for a continuous-time system `sys`, else returns `false`."""
iscontinuous(sys::LTISystem) = timeevol(sys) isa Continuous
"""`isdiscrete(sys)`

Returns `true` for a discrete-time system `sys`, else returns `false`."""
isdiscrete(sys::LTISystem) = timeevol(sys) isa Discrete


function Base.getproperty(sys::LTISystem, s::Symbol)
    if s === :Ts
        # if !isdiscrete(sys) # NOTE this line seems to be breaking inference of isdiscrete (is there a test for this?)
        if isdiscrete(sys)
            return timeevol(sys).Ts
        else
            @warn "Getting time 0.0 for non-discrete systems is deprecated. Check `isdiscrete` before trying to access time."
            return 0.0
        end
    else
        return getfield(sys, s)
    end
end

Base.propertynames(sys::LTISystem, private::Bool=false) =
    (fieldnames(typeof(sys))..., :nu, :ny, (isdiscrete(sys) ? (:Ts,) : ())...)




common_timeevol(systems::LTISystem...) = common_timeevol(timeevol(sys) for sys in systems)


"""
    isstable(sys)

Returns `true` if `sys` is stable, else returns `false`."""
isstable(sys::LTISystem{Continuous}) = all(real.(poles(sys)) .< 0)
isstable(sys::LTISystem{<:Discrete}) = all(abs.(poles(sys)) .< 1)

# Fallback since LTISystem not AbstractArray
Base.size(sys::LTISystem, i::Integer) = size(sys)[i]


## Names =======================================================================

"""
    system_name(nothing::LTISystem)

Return the name of the system. If the system does not have a name, an empty string is returned.
"""
system_name(::LTISystem) = ""

"""
    input_names(P)
    input_names(P, i)

Get a vector of strings with the names of the inputs of `P`, or the `i`:th name if an index is given.
"""
input_names(P::LTISystem; kwargs...) = [input_names(P, i; kwargs...) for i in 1:ninputs(P)]
function input_names(P::LTISystem, i; default = "u")
    i <= ninputs(P) || throw(BoundsError(P, i))
    ninputs(P) == 1 && (return default)
    "u($i)"
end

"""
    output_names(P)
    output_names(P, i)

Get a vector of strings with the names of the outputs of `P`, or the `i`:th name if an index is given.
"""
output_names(P::LTISystem; kwargs...) = [output_names(P, i; kwargs...) for i in 1:noutputs(P)]
function output_names(P::LTISystem, i; default = "y")
    i <= noutputs(P) || throw(BoundsError(P, i))
    noutputs(P) == 1 && (return default)
    "y($i)"
end

"""
    state_names(P)
    state_names(P, i)

Get a vector of strings with the names of the states of `P`, or the `i`:th name if an index is given.
"""
state_names(P::LTISystem; kwargs...) = [state_names(P, i; kwargs...) for i in 1:nstates(P)]
function state_names(P::LTISystem, i; default = "x")
    i <= nstates(P) || throw(BoundsError(P, i))
    nstates(P) == 1 && (return default)
    "x($i)"
end