struct NamedSystem{T<:LTISystem} <: LTISystem
    sys::T
    inputnames::Vector{String}
    outputnames::Vector{String}
    function NamedSystem{T}(sys::T, inputnames::Vector{String}, outputnames::Vector{String}) where {T<:LTISystem}
        ny,nu = size(sys)
        @assert nu == length(inputnames)             "Must have same number of inputnames as inputs"
        @assert ny == length(outputnames)            "Must have same number of outputnames as outputs"
        new{T}(sys, inputnames, outputnames)
    end
end

### Fallbacks for LTI system
inputnames(sys::LTISystem)  = fill("", size(sys,2))
outputnames(sys::LTISystem) = fill("", size(sys,1))

### For NamedSystem
inputnames(sys::NamedSystem)  = sys.inputnames
outputnames(sys::NamedSystem) = sys.outputnames

function NamedSystem(sys::T; inputnames="", outputnames="") where {T}
    inputnames  = validate_names(inputnames,  "inputnames",  size(sys,2))
    outputnames = validate_names(outputnames, "outputnames", size(sys,1))
    NamedSystem{T}(sys, inputnames, outputnames)
end

function validate_names(names, nametype, n)
    if names == ""
        return fill("", n)
    elseif isa(names, Vector) && eltype(names) <: AbstractString
        return String[names[i] for i = 1:n]
    elseif isa(names, AbstractString)
        return String[names * "$i" for i = 1:n]
    else
        error("$nametype must be of type `AbstractString` or Vector{AbstractString}")
    end
end

# Used in @show
function format_names(names::Vector{String}, default::AbstractString, unknown::AbstractString)
    n = size(names, 1)
    if all(names .== "")
        return String[default * string(i) for i=1:n]
    else
        for (i, name) in enumerate(names)
            names[i] = String((name == "") ? unknown : n)
        end
        return names
    end
end

function append(systems::NamedSystem...)
     outputs = vcat([s.outputnames for s in systems]...)
     inputs = vcat([s.inputnames for s in systems]...)
    NamedSystem(append(systems...), inputs, output)
end

function ==(s1::NamedSystem, s2::NamedSystem)
    if s1.sys == s2.sys
        fields = [:inputnames, :outputnames]
        for field in fields
            if getfield(s1, field) != getfield(s2, field)
                return false
            end
        end
        return true
    else
        return false
    end
end

function +(s1::NamedSystem, s2::NamedSystem)
    sys = s1.sys + s2.sys
    # Naming strategy: If only one sys is named, use that. If the names are the
    # same, use them. If the names conflict, then they are ignored, and the
    # default "" is used.
    if all(s1.inputnames .== "")
        inputnames = s2.inputnames
    elseif all(s2.inputnames .== "") || (s1.inputnames == s2.inputnames)
        inputnames = s1.inputnames
    else
        inputnames = fill(String(""), size(s1,1))
    end
    if all(s1.outputnames .== "")
        outputnames = s2.outputnames
    elseif all(s2.outputnames .== "") || (s1.outputnames == s2.outputnames)
        outputnames = s1.outputnames
    else
        outputnames = fill(String(""), size(s1,2))
    end
    return NamedSystem(sys, inputnames, outputnames)
end

+(s::NamedSystem, n::Real) = NamedSystem(s.sys+n, s.inputnames, s.outputnames)

-(s::NamedSystem) = NamedSystem(-s.sys, s.inputnames, s.outputnames)

*(s1::NamedSystem, s2::NamedSystem) = NamedSystem(s1.sys*s2.sys, s2.inputnames, s1.outputnames)

/(n::Real, s::NamedSystem) = NamedSystem(n/s.sys, s.outputnames, s.inputnames)

/(s::NamedSystem, n::Real) = NamedSystem(s.sys/n, s.inputnames, s.outputnames)

Base.getindex(s::NamedSystem, inds...) =
    NamedSystem(getindex(s, inds...), [s.inputnames[cols];], [s.outputnames[rows];])
