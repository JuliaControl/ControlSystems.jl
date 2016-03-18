# Model interconnections

series(s1::LTISystem, s2::LTISystem) = s2*s1

parallel(s1::LTISystem, s2::LTISystem) = s1 + s2

append() = LTISystem[]
function append(systems::StateSpace...)
    Ts = systems[1].Ts
    if !all([s.Ts == Ts for s in systems])
        error("Sampling time mismatch")
    end
    A = blkdiag([s.A for s in systems]...)
    B = blkdiag([s.B for s in systems]...)
    C = blkdiag([s.C for s in systems]...)
    D = blkdiag([s.D for s in systems]...)
    states = vcat([s.statenames for s in systems]...)
    outputs = vcat([s.outputnames for s in systems]...)
    inputs = vcat([s.inputnames for s in systems]...)
    return StateSpace(A, B, C, D, Ts, states, inputs, outputs)
end

function append(systems::TransferFunction...)
    Ts = systems[1].Ts
    if !all([s.Ts == Ts for s in systems])
        error("Sampling time mismatch")
    end
    mat = blkdiag([s.matrix for s in systems]...)
    inputs = vcat([s.inputnames for s in systems]...)
    outputs = vcat([s.outputnames for s in systems]...)
    return TransferFunction(mat, Ts, inputs, outputs)
end

append(systems::LTISystem...) = append(promote(systems...)...)

#This is needed until julia removes deprecated vect()
function Base.vect(sys::Union{LTISystem, Real}...)
    T = Base.promote_typeof(sys...)
    copy!(Array(T,length(sys)), sys)
end

function Base.vcat(systems::StateSpace...)
    # Perform checks
    nu = systems[1].nu
    if !all([s.nu == nu for s in systems])
        error("All systems must have same input dimension")
    end
    Ts = systems[1].Ts
    if !all([s.Ts == Ts for s in systems])
        error("Sampling time mismatch")
    end
    A = blkdiag([s.A for s in systems]...)
    B = vcat([s.B for s in systems]...)
    C = blkdiag([s.C for s in systems]...)
    D = vcat([s.D for s in systems]...)
    states = vcat([s.statenames for s in systems]...)
    outputs = vcat([s.outputnames for s in systems]...)
    inputs = systems[1].inputnames
    if !all([s.inputnames == inputs for s in systems])
        inputs = fill(UTF8String(""),size(inputs, 1))
    end
    return StateSpace(A, B, C, D, Ts, states, inputs, outputs)
end

function Base.vcat(systems::TransferFunction...)
    # Perform checks
    nu = systems[1].nu
    if !all([s.nu == nu for s in systems])
        error("All systems must have same input dimension")
    end
    Ts = systems[1].Ts
    if !all([s.Ts == Ts for s in systems])
        error("Sampling time mismatch")
    end
    mat = vcat([s.matrix for s in systems]...)
    outputs = vcat([s.outputnames for s in systems]...)
    inputs = systems[1].inputnames
    if !all([s.inputnames == inputs for s in systems])
        inputs = fill(UTF8String(""),size(inputs, 1))
    end
    return TransferFunction(mat, Ts, inputs, outputs)
end

Base.vcat(systems::LTISystem...) = vcat(promote(systems...)...)

function Base.vcat{T<:Real}(systems::Union{VecOrMat{T},T,TransferFunction}...)
    if promote_type(map(e->typeof(e),systems)...) <: TransferFunction
        vcat(map(e->convert(TransferFunction,e),systems)...)
    else
        cat(1,systems...)
    end
end

function Base.hcat(systems::StateSpace...)
    # Perform checks
    ny = systems[1].ny
    if !all([s.ny == ny for s in systems])
        error("All systems must have same output dimension")
    end
    Ts = systems[1].Ts
    if !all([s.Ts == Ts for s in systems])
        error("Sampling time mismatch")
    end
    A = blkdiag([s.A for s in systems]...)
    B = blkdiag([s.B for s in systems]...)
    C = hcat([s.C for s in systems]...)
    D = hcat([s.D for s in systems]...)
    states = vcat([s.statenames for s in systems]...)
    inputs = vcat([s.inputnames for s in systems]...)
    outputs = systems[1].outputnames
    if !all([s.outputnames == outputs for s in systems])
        outputs = fill(UTF8String(""),size(outputs, 1))
    end
    return StateSpace(A, B, C, D, Ts, states, inputs, outputs)
end

function Base.hcat(systems::TransferFunction...)
    # Perform checks
    ny = systems[1].ny
    if !all([s.ny == ny for s in systems])
        error("All systems must have same output dimension")
    end
    Ts = systems[1].Ts
    if !all([s.Ts == Ts for s in systems])
        error("Sampling time mismatch")
    end
    mat = hcat([s.matrix for s in systems]...)
    inputs = vcat([s.inputnames for s in systems]...)
    outputs = systems[1].outputnames
    if !all([s.outputnames == outputs for s in systems])
        outputs = fill(UTF8String(""),size(outputs, 1))
    end
    return TransferFunction(mat, Ts, inputs, outputs)
end

Base.hcat(systems::LTISystem...) = hcat(promote(systems...)...)

function Base.hcat{T<:Real}(systems::Union{T,VecOrMat{T},TransferFunction}...)
    if promote_type(map(e->typeof(e),systems)...) <: TransferFunction
        hcat(map(e->convert(TransferFunction,e),systems)...)
    else
        cat(2,systems...)
    end
end

# Empty definition to get rid of warning
Base.blkdiag() = []
function Base.blkdiag(mats::Matrix...)
    rows = Int[size(m, 1) for m in mats]
    cols = Int[size(m, 2) for m in mats]
    T = eltype(mats[1])
    for ind=1:length(mats)
        T = promote_type(T, eltype(mats[ind]))
    end
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
