
struct TimedVal{ValT, TimeT}
    t::TimeT
    val::ValT
end

# TODO: Create index of recently accessed times for quicker search
struct DDEBuffer{ValT, TimeT}
    buff::Vector{Deque{TimedVal{ValT,TimeT}}}
    τ::Vector{TimeT}
    function DDEBuffer(::Type{ValT}, τ::AbstractVector{TimeT}) where {ValT,TimeT}
        buff = [Deque{TimedVal{ValT,TimeT}}() for i in 1:length(τ)]
        new{ValT,TimeT}(buff, copy(τ))
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
    push!(b.buff[indx], TimedVal(t, val))
    
    # println(length(b.buff[indx]))
end

# Get last value before time t
function (b::DDEBuffer)(t, indx::Int)
    bi = b.buff[indx]
    val = first(bi).val
    if first(bi).t > t
        error("time $t is before start of buffer")
    end
    for vt in bi
        if vt.t <= t
            val = vt.val
        else
            break
        end
    end
    return val
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

function Base.setindex!(hf::HistoryFunction, val, t, indx::Int)
    # Push to buffer
    push!(hf.buff, t, val, indx)
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


