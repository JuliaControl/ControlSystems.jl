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
        inputs = UTF8String["" for i = 1:size(inputs, 1)]
    end
    return StateSpace(A, B, C, D, Ts, states, inputs, outputs)
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
        outputs = UTF8String["" for i = 1:size(outputs, 1)]
    end
    return StateSpace(A, B, C, D, Ts, states, inputs, outputs)
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
        inputs = UTF8String["" for i = 1:size(inputs, 1)]
    end
    return TransferFunction(mat, Ts, inputs, outputs)
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
        outputs = UTF8String["" for i = 1:size(outputs, 1)]
    end
    return TransferFunction(mat, Ts, inputs, outputs)
end
