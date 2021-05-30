abstract type AbstractResult end # Common for all result types, e.g., SimResult and FreqRespResult

## SimResult: the output of lsim etc. ==========================================
abstract type AbstractSimResult <: AbstractResult end # Result of a time-domain simulation

struct SimResult{Ty, Tt, Tx, Tu, Ts} <: AbstractSimResult # Result of lsim
    y::Ty
    t::Tt
    x::Tx
    u::Tu
    sys::Ts
end

# To emulate, e.g., lsim(sys, u)[1] -> y
function Base.getindex(r::SimResult, i::Int)
    return getfield(r, i) 
end

function Base.getindex(r::SimResult, v::AbstractVector)
    return getfield.((r,), v) 
end

# to allow destructuring, e.g., y,t,x = lsim(sys, u)
# This performs explicit iteration in the type domain to ensure inferability
Base.iterate(r::SimResult)              = (r.y, Val(:t))
Base.iterate(r::SimResult, ::Val{:t})   = (r.t, Val(:x))
Base.iterate(r::SimResult, ::Val{:x})   = (r.x, Val(:u))
Base.iterate(r::SimResult, ::Val{:u})   = (r.u, Val(:done))
Base.iterate(r::SimResult, ::Val{:done}) = nothing


function Base.getproperty(r::SimResult, s::Symbol)
    s ∈ fieldnames(SimResult) && (return getfield(r, s))
    if s === :nx
        return size(r.x, 1)
    elseif s === :nu
        return size(r.u, 1)
    elseif s === :ny
        return size(r.y, 1)
    else
        throw(ArgumentError("Unsupported property $s"))
    end
end