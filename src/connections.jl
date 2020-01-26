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
function append(systems::(ST where ST<:AbstractStateSpace)...)
    ST = promote_type(typeof.(systems)...)
    Ts = common_sample_time(s.Ts for s in systems)
    A = blockdiag([s.A for s in systems]...)
    B = blockdiag([s.B for s in systems]...)
    C = blockdiag([s.C for s in systems]...)
    D = blockdiag([s.D for s in systems]...)
    return ST(A, B, C, D, Ts)
end

function append(systems::TransferFunction...)
    Ts = common_sample_time(s.Ts for s in systems)
    mat = blockdiag([s.matrix for s in systems]...)
    return TransferFunction(mat, Ts)
end

append(systems::LTISystem...) = append(promote(systems...)...)


function Base.vcat(systems::DelayLtiSystem...)
    P = vcat_1([sys.P for sys in systems]...) # See PartitionedStateSpace
    Tau = vcat([sys.Tau for sys in systems]...)
    return DelayLtiSystem(P, Tau)
end

function Base.hcat(systems::DelayLtiSystem...)
    P = hcat_1([sys.P for sys in systems]...)  # See PartitionedStateSpace
    Tau = vcat([sys.Tau for sys in systems]...)
    return DelayLtiSystem(P, Tau)
end


function Base.vcat(systems::ST...) where ST <: AbstractStateSpace
    # Perform checks
    nu = systems[1].nu
    if !all(s.nu == nu for s in systems)
        error("All systems must have same input dimension")
    end
    A = blockdiag([s.A for s in systems]...)
    B = vcat([s.B for s in systems]...)
    C = blockdiag([s.C for s in systems]...)
    D = vcat([s.D for s in systems]...)
    Ts = common_sample_time(s.Ts for s in systems)
    return ST(A, B, C, D, Ts)
end

function Base.vcat(systems::TransferFunction...)
    # Perform checks
    nu = systems[1].nu
    if !all(s.nu == nu for s in systems)
        error("All systems must have same input dimension")
    end
    Ts = common_sample_time(s.Ts for s in systems)
    mat = vcat([s.matrix for s in systems]...)

    return TransferFunction(mat, Ts)
end

Base.vcat(systems::LTISystem...) = vcat(promote(systems...)...)

function Base.hcat(systems::ST...) where ST <: AbstractStateSpace
    # Perform checks
    ny = systems[1].ny
    if !all(s.ny == ny for s in systems)
        error("All systems must have same output dimension")
    end
    Ts = common_sample_time(s.Ts for s in systems)
    A = blockdiag([s.A for s in systems]...)
    B = blockdiag([s.B for s in systems]...)
    C = hcat([s.C for s in systems]...)
    D = hcat([s.D for s in systems]...)

    return ST(A, B, C, D, Ts)
end

function Base.hcat(systems::TransferFunction...)
    # Perform checks
    ny = systems[1].ny
    if !all(s.ny == ny for s in systems)
        error("All systems must have same output dimension")
    end
    Ts = common_sample_time(s.Ts for s in systems)
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
function Base.typed_hcat(::Type{T}, X...) where {T<:LTISystem}
    hcat(convert.(T, X)...)
end
# Ambiguity
Base.typed_hcat(::Type{S}, X::Number...) where {S<:LTISystem} = hcat(convert.(S, X)...)
Base.typed_hcat(::Type{S}, X::Union{AbstractArray{<:Number,1}, AbstractArray{<:Number,2}}...) where {S<:LTISystem} = hcat(convert.(S, X)...)

# Catch special cases where inv(sys) might not be possible after promotion, like improper tf
function /(sys1::Union{StateSpace,AbstractStateSpace}, sys2::LTISystem)
    sys1new, sys2new = promote(sys1, 1/sys2)
    return sys1new*sys2new
end

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
