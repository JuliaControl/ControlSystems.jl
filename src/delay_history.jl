struct TimedVal{ValT, TimeT}
    t::TimeT
    val::ValT
end

# TODO: Create index of recently accessed times for quicker search
struct DDEBuffer{ValT, TimeT}
    buff::Vector{Deque{TimedVal{ValT,TimeT}}}
    τ::Vector{TimeT}
    disc::Vector{ValT}
    function DDEBuffer(::Type{ValT}, τ::AbstractVector{TimeT}) where {ValT,TimeT}
        buff = [Deque{TimedVal{ValT,TimeT}}() for i in 1:length(τ)]
        new{ValT,TimeT}(buff, copy(τ), ValT[])
    end
end

# get the element at the back
Base.last(b::DDEBuffer, indx) = last(b.buff[indx])
# test whether the buffer is empty
Base.isempty(b::DDEBuffer, indx::Int) = isempty(b.buff[indx])
# get the number of elements
Base.length(b::DDEBuffer, indx::Int) = length(b.buff[indx])  
# remove an element from the back
Base.pop!(b::DDEBuffer, indx::Int) = pop!(b.buff[indx])  
# add an element to the front
Base.pushfirst!(b::DDEBuffer, val, indx::Int) = pushfirst!(b.buff[indx], val)
# remove an element from the front
Base.popfirst!(b::DDEBuffer, indx::Int) = popfirst!(b.buff[indx])  
# get the element at the front
Base.first(b::DDEBuffer, indx::Int) = first(b.buff[indx])  

# add an element to the back
function Base.push!(b::DDEBuffer, t, val, indx::Int)
    if !isempty(b,indx) && last(b, indx).t > t
        throw(AssertionError("Can only push times after the buffer end"))
    end
    push!(b.buff[indx], TimedVal(t, val))
end

# Get last value before time t
function (b::DDEBuffer)(t, indx::Int)
    bi = b.buff[indx]
    v0 = nothing
    t0 = nothing
    v1 = nothing
    t1 = nothing
    v2 = nothing
    t2 = nothing

    if first(bi).t > t
        throw(AssertionError("time $t is before start of buffer"))
    end

    for vt in bi
        if vt.t <= t
            # Keep updating as long as vt.t < t
            # Save prev value
            v0 = v1
            t0 = t1
            # Current value
            v1 = vt.val
            t1 = vt.t
        else
            # We found te last value, now get one more
            v2 = vt.val
            t2 = vt.t
            break

        end
    end
    # Error if no value found
    v1 === nothing && throw(AssertionError("Value not in buffer range at t: $t"))

    if t1 == t
        # If exact value
        return v1
    elseif v2 === nothing
        # If extrapolating forward
        @warn "Interpolation inaccurate at t=$t"
        return v1
    elseif t2 <= t1
        @warn "Interpolation found two equal time steps at t=$t"
        return v1
    elseif t2 in b.disc
        # Next time point is a discontinuity
        if v0 === nothing
            # Not no provious point exists
            @warn "Discontinuity found at next time point, no interpolation t=$t, t2=$t2"
            return v1
        elseif t0 == t1
            # Previous point is same as current
            @warn "Discontinuity found at next time point, two equal vals t=$t, t2=$t2"
            return v1
        else
            # Extrapolate linearly from two previous points
            return (v0*(t1 - t) + v1*(t - t0))/(t1-t0)
        end
    else
        # Interpolate between current and next
        return (v1*(t2 - t) + v2*(t - t1))/(t2-t1)
    end
    return v1 # Can't happen
end

"""

"""
struct HistoryFunction{ValT,TimeT,TF<:Function}
    buff::DDEBuffer{ValT, TimeT}
    tinit::TimeT
    history::TF
    eps::TimeT
end

"""
    `h = HistoryFunction(::Type{ValT}, f::Function, tinit, τ; eps=sqrt(ValT)) where ValT`

    Returns object `h` that is able to store and return a history of old values.
    It will store up to `length(τ)` values.

    h(t,i) will return the stored value with index `i` at a previous time `t`.
    If `t < tinit` then h(t,i) = f(t,i).

    `h[t,i] = v` Saves the value `v` for index `i` at time `t.
    It will also erase old values at index `i` for times before `t-τ[i]` 

    `eps` is the extra time old values will be saved.
"""
function HistoryFunction(::Type{ValT}, f::Function, tinit, τ; eps=sqrt(eps(ValT))) where ValT
    buff = DDEBuffer(ValT, τ)
    HistoryFunction(buff, tinit, f, eps)
end

function (hf::HistoryFunction)(t, indx::Int)
    if t < hf.tinit
        return hf.history(t, indx)
    else
        return hf.buff(t, indx)
    end
end

function Base.setindex!(hf::HistoryFunction, val, t, indx::Int, integrator)
    # Push to buffer
    push!(hf.buff, t, val, indx)
    # Add potential discontinuity to list
    if !isempty(integrator.opts.d_discontinuities) && first(integrator.opts.d_discontinuities) == t
        !(t in hf.buff.disc) && push!(hf.buff.disc, t)
    end
    # Erase history before t-τ
    t0 = t-hf.buff.τ[indx] - hf.eps
    
    if first(hf.buff, indx).t < t0
        # Start erasing
        while true
            v = popfirst!(hf.buff, indx)
            # Make sure second value is still in range, since we search forward
            if first(hf.buff, indx).t < t0
                continue # erasing
            else
                # Push it back
                pushfirst!(hf.buff, v, indx)
                break # no more to erase
            end
        end
    end
end


