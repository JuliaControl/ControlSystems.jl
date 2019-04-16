# Model interconnections

"""
`series(sys1::LTISystem, sys2::LTISystem)`

Connect systems in series, equivalent to `sys2*sys1`
"""
series(sys1::LTISystem, sys2::LTISystem) = sys2*sys1

"""
`series(sys1::LTISystem, sys2::LTISystem)`

Connect systems in parallel, equivalent to `sys2+sys1`
"""
parallel(sys1::LTISystem, sys2::LTISystem) = sys1 + sys2

append() = LTISystem[]
"""
`append(systems::StateSpace...), append(systems::TransferFunction...)`

Append systems in block diagonal form
"""
function append(systems::StateSpace...)
    Ts = systems[1].Ts
    if !all(s.Ts == Ts for s in systems)
        error("Sampling time mismatch")
    end
    A = blockdiag([s.A for s in systems]...)
    B = blockdiag([s.B for s in systems]...)
    C = blockdiag([s.C for s in systems]...)
    D = blockdiag([s.D for s in systems]...)
    return StateSpace(A, B, C, D, Ts)
end

function append(systems::TransferFunction...)
    Ts = systems[1].Ts
    if !all(s.Ts == Ts for s in systems)
        error("Sampling time mismatch")
    end
    mat = blockdiag([s.matrix for s in systems]...)
    return TransferFunction(mat, Ts)
end

append(systems::LTISystem...) = append(promote(systems...)...)

function Base.vcat(systems::StateSpace...)
    # Perform checks
    nu = systems[1].nu
    if !all(s.nu == nu for s in systems)
        error("All systems must have same input dimension")
    end
    Ts = systems[1].Ts
    if !all(s.Ts == Ts for s in systems)
        error("Sampling time mismatch")
    end
    A = blockdiag([s.A for s in systems]...)
    B = vcat([s.B for s in systems]...)
    C = blockdiag([s.C for s in systems]...)
    D = vcat([s.D for s in systems]...)

    return StateSpace(A, B, C, D, Ts)
end

function Base.vcat(systems::TransferFunction...)
    # Perform checks
    nu = systems[1].nu
    if !all(s.nu == nu for s in systems)
        error("All systems must have same input dimension")
    end
    Ts = systems[1].Ts
    if !all(s.Ts == Ts for s in systems)
        error("Sampling time mismatch")
    end
    mat = vcat([s.matrix for s in systems]...)
    return TransferFunction(mat, Ts)
end

Base.vcat(systems::LTISystem...) = vcat(promote(systems...)...)

function Base.hcat(systems::StateSpace...)
    # Perform checks
    ny = systems[1].ny
    if !all(s.ny == ny for s in systems)
        error("All systems must have same output dimension")
    end
    Ts = systems[1].Ts
    if !all(s.Ts == Ts for s in systems)
        error("Sampling time mismatch")
    end
    A = blockdiag([s.A for s in systems]...)
    B = blockdiag([s.B for s in systems]...)
    C = hcat([s.C for s in systems]...)
    D = hcat([s.D for s in systems]...)

    return StateSpace(A, B, C, D, Ts)
end

function Base.hcat(systems::TransferFunction...)
    # Perform checks
    ny = systems[1].ny
    if !all(s.ny == ny for s in systems)
        error("All systems must have same output dimension")
    end
    Ts = systems[1].Ts
    if !all(s.Ts == Ts for s in systems)
        error("Sampling time mismatch")
    end
    mat = hcat([s.matrix for s in systems]...)
    return TransferFunction(mat, Ts)
end

Base.hcat(systems::LTISystem...) = hcat(promote(systems...)...)

function Base._cat_t(::Val{1}, T::Type{<:LTISystem}, X...)
        vcat(convert.(T, X)...)
end

function Base._cat_t(::Val{2}, T::Type{<:LTISystem}, X...)
        hcat(convert.(T, X)...)
end

# Used in typed_hvcat
Base.typed_hcat(::Type{T}, X...) where {T<:LTISystem} = hcat(convert.(T, X)...)
# Ambiguity
Base.typed_hcat(::Type{T}, X::Number...) where {T<:LTISystem, N} = hcat(convert.(T, X)...)

# function hvcat(rows::Tuple{Vararg{Int}}, systems::Union{Number,AbstractVecOrMat{<:Number},LTISystem}...)
#     T = Base.promote_typeof(systems...)
#     nbr = length(rows)  # number of block rows
#     rs = Array{T,1}(nbr)
#     a = 1
#     for i = 1:nbr
#         rs[i] = hcat(convert.(T,systems[a:a-1+rows[i]])...)
#         a += rows[i]
#     end
#     vcat(rs...)
# end

# function _get_common_sampling_time(sys_vec::Union{AbstractVector{LTISystem},AbstractVecOrMat{<:Number},Number})
#     Ts = -1.0 # Initalize corresponding to undefined sampling time
#
#     for sys in sys_vec
#         if !all(s.Ts == Ts for s in systems])
#             error("Sampling time mismatch")
#         end
#     end
#
# end


# function Base.hcat{T<:Number}(systems::Union{T,AbstractVecOrMat{T},TransferFunction}...)
#     S = promote_type(map(e->typeof(e),systems)...) # TODO: Should be simplified
#
#     idx_first_tf = findfirst(e -> isa(e, TransferFunction), systems)
#     Ts = sys_tuple[idx_first_tf].Ts
#
#     if S <: TransferFunction
#         hcat(map(e->convert(TransferFunction,e),systems)...)
#     else
#         cat(2,systems...)
#     end
# end


blockdiag(mats::AbstractMatrix...) = blockdiag(promote(mats...)...)

function blockdiag(mats::AbstractMatrix{T}...) where T
    rows = Int[size(m, 1) for m in mats]
    cols = Int[size(m, 2) for m in mats]
    res = zeros(T, sum(rows), sum(cols))
    m = 1
    n = 1
    for ind=1:length(mats)
        mat = mats[ind]
        i = rows[ind]
        j = cols[ind]
        res[m:m + i - 1, n:n + j - 1] = mat
        m += i
        n += j
    end
    return res
end
