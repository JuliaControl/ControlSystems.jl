abstract type AbstractResult end # Common for all result types, e.g., SimResult and FreqRespResult

## SimResult: the output of lsim etc. ==========================================
abstract type AbstractSimResult <: AbstractResult end # Result of a time-domain simulation


"""
    SimResult{Ty, Tt, Tx, Tu, Ts} <: AbstractSimResult

Result structure containing the results of time-domain simulations using `lsim, step, impulse`.
The structure can be plotted using
```julia
result = lsim(...)
plot(result, plotu=true, plotx=false)
```
and destructured like
```julia
y, t, x, u = result
```

# Fields:
- `y::Ty`
- `t::Tt`
- `x::Tx`
- `u::Tu`
- `sys::Ts`

## Concatenation of SimResults

Two SimResults can be concatenated in time using `hcat`, or `[res1 res2]`, the rules for this are as follows:
- If the start time of the second result is one sample interval after the end time of the first result, the results are concatenated and the length of the result is the sum of the lengths of the two results.
- If the start time of the second result is equal to the end time of the first result, _and_ the initial state of the second result is equal to the final state of the first result, the results are concatenated omitting the initial point from the second result, which would otherwise have been repeated. The length of the result is the sum of the lengths of the two results minus one.
If none of the above holds, a warning is issued and the result has the length of the sum of the lengths of the two results.
- If the sample intervals of the two results are different, an error is thrown.
"""
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

function Base.hcat(r1::SimResult, r2::SimResult)
    r1.sys.Ts == r2.sys.Ts || throw(ArgumentError("Sampling-time mismatch"))
    if r1.x[:, end] == r2.x[:, 1] && r1.t[end] == r2.t[1]
        r1.u[:, end] == r2.u[:, 1] || @warn "Concatenated SimResults have different inputs at the join"
        # This is a common case when r2 starts with initial conditions from r1
        return SimResult(hcat(r1.y, r2.y[:, 2:end]), vcat(r1.t, r2.t[2:end]), hcat(r1.x, r2.x[:, 2:end]), hcat(r1.u, r2.u[:, 2:end]), r1.sys)
    elseif !(r1.t[end] + r1.sys.Ts ≈ r2.t[1])
        @warn "Concatenated SimResults do not appear to be continuous in time, the first ends at t=$(r1.t[end]) and the second starts at t=$(r2.t[1]). With sample interval Ts=$(r1.sys.Ts), the second simulation was expected to start at t=$(r1.t[end] + r1.sys.Ts) To start a simulation at a non-zero time, pass a time vector to lsim."
    end
    SimResult(hcat(r1.y, r2.y), vcat(r1.t, r2.t), hcat(r1.x, r2.x), hcat(r1.u, r2.u), r1.sys)
end

issiso(r::SimResult) = issiso(r.sys)

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


"""
    RootLocusResult{T} <: AbstractResult

Result structure containing the result of the root locus using `rlocus`.
The structure can be plotted using
```julia
result = rlocus(...)
plot(result; plotu=true, plotx=false)?
```
and destructured like
```julia
roots, Z, K = result
```

# Fields:
- `roots::Tr`
- `Z::Tz`
- `K::Tk`
- `sys::Ts`
"""
struct RootLocusResult{Tr, Tz, Tk, Ts} <: AbstractResult # Result of rlocus
    roots::Tr
    Z::Tz
    K::Tk
    sys::Ts
end

Base.iterate(r::RootLocusResult)               = (r.roots, Val(:Z))
Base.iterate(r::RootLocusResult, ::Val{:Z})    = (r.Z, Val(:K))
Base.iterate(r::RootLocusResult, ::Val{:K})    = (r.K, Val(:done))
Base.iterate(r::RootLocusResult, ::Val{:done}) = nothing
