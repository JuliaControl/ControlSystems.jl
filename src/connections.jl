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
    Ts = systems[1].Ts
    if !all(s.Ts == Ts for s in systems)
        error("Sampling time mismatch")
    end
    A = blockdiag([s.A for s in systems]...)
    B = blockdiag([s.B for s in systems]...)
    C = blockdiag([s.C for s in systems]...)
    D = blockdiag([s.D for s in systems]...)
    return ST(A, B, C, D, Ts)
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
    Ts = systems[1].Ts
    if !all(s.Ts == Ts for s in systems)
        error("Sampling time mismatch")
    end
    A = blockdiag([s.A for s in systems]...)
    B = vcat([s.B for s in systems]...)
    C = blockdiag([s.C for s in systems]...)
    D = vcat([s.D for s in systems]...)

    return ST(A, B, C, D, Ts)
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

function Base.hcat(systems::ST...) where ST <: AbstractStateSpace
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

    return ST(A, B, C, D, Ts)
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
function Base.typed_hcat(::Type{T}, X...) where {T<:LTISystem}
    hcat(convert.(T, X)...)
end
# Ambiguity
Base.typed_hcat(::Type{T}, X::Number...) where {T<:LTISystem, N} = hcat(convert.(T, X)...)

# Catch special cases where inv(sys) might not be possible after promotion, like improper tf
function /(sys1::Union{StateSpace,AbstractStateSpace}, sys2::LTISystem)
    sys1new, sys2new = promote(sys1, 1/sys2)
    return sys1new*sys2new
end

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



"""
`feedback(L)` Returns L/(1+L)
`feedback(P1,P2)` Returns P1/(1+P1*P2)
"""
feedback(L::TransferFunction) = L/(1+L)
feedback(P1::TransferFunction, P2::TransferFunction) = P1/(1+P1*P2)

#Efficient implementations
function feedback(L::TransferFunction{T}) where T<:SisoRational
    if size(L) != (1,1)
        error("MIMO TransferFunction feedback isn't implemented, use L/(1+L)")
    end
    P = numpoly(L)
    Q = denpoly(L)
    tf(P, P+Q, L.Ts)
end

function feedback(L::TransferFunction{T}) where {T<:SisoZpk}
    if size(L) != (1,1)
        error("MIMO TransferFunction feedback isn't implemented, use L/(1+L)")
    end
    #Extract polynomials and create P/(P+Q)
    k = L.matrix[1].k
    denpol = numpoly(L)[1]+denpoly(L)[1]
    kden = denpol[end] # Get coeff of s^n
    # Create siso system
    sisozpk = T(L.matrix[1].z, roots(denpol), k/kden)
    return TransferFunction{T}(fill(sisozpk,1,1), L.Ts)
end

"""
`feedback(sys)`

`feedback(sys1,sys2)`

Forms the negative feedback interconnection
```julia
>-+ sys1 +-->
  |      |
 (-)sys2 +
```
If no second system is given, negative identity feedback is assumed
"""
function feedback(sys::Union{StateSpace, DelayLtiSystem})
    ninputs(sys) != noutputs(sys) && error("Use feedback(sys1, sys2) if number of inputs != outputs")
    feedback(sys,ss(Matrix{numeric_type(sys)}(I,size(sys)...)))
end

function feedback(sys1::StateSpace,sys2::StateSpace)
    sum(abs.(sys1.D)) != 0 && sum(abs.(sys2.D)) != 0 && error("There can not be a direct term (D) in both sys1 and sys2")
    A = [sys1.A+sys1.B*(-sys2.D)*sys1.C sys1.B*(-sys2.C); sys2.B*sys1.C  sys2.A+sys2.B*sys1.D*(-sys2.C)]
    B = [sys1.B; sys2.B*sys1.D]
    C = [sys1.C  sys1.D*(-sys2.C)]
    ss(A,B,C,sys1.D)
end


"""
    feedback(s1::AbstractStateSpace, s2::AbstractStateSpace;
        U1=1:size(s1,2), Y1=1:size(s1,1), U2=1:size(s2,2), Y2=1:size(s2,1),
        W1=1:size(s1,2), Z1=1:size(s1,1), W2=Int[], Z2=Int[],
        w_reorder=nothing, z_reorder=nothing,
        neg_feeback::Bool=true)


`U1`, `Y1`, `U2`, `Y2` contains the indices of the signals that should be connected.
`W1`, `Z1`, `W2`, `Z2` contains the signal indices of `s1` and `s2` that should be kept.

Specify  `Wreorder` and `Zreorder` to reorder [W1; W2] and [Z1; Z2] after matrix multiplications.

Negative feedback is the default. Specify `neg_feedback=false` for positive feedback.

See Zhou etc. for similar (somewhat less symmetric) formulas.
"""
function feedback(s1::AbstractStateSpace, s2::AbstractStateSpace;
    U1=1:size(s1,2), Y1=1:size(s1,1), U2=1:size(s2,2), Y2=1:size(s2,1),
    W1=1:size(s1,2), Z1=1:size(s1,1), W2=Int[], Z2=Int[],
    w_reorder=nothing, z_reorder=nothing,
    neg_feedback::Bool=true)

    if !(allunique(Y1) && allunique(Y2))
        warning("Connecting a single output to multiple inputs.")
    end
    if !(allunique(U1) && allunique(U2))
        warning("Connecting multiple outputs to a single input.")
    end

    α = neg_feedback ? -1 : 1 # The sign of feedback

    s1_B1 = s1.B[:,W1]
    s1_B2 = s1.B[:,U1]
    s1_C1 = s1.C[Z1,:]
    s1_C2 = s1.C[Y1,:]
    s1_D11 = s1.D[Z1,W1]
    s1_D12 = s1.D[Z1,U1]
    s1_D21 = s1.D[Y1,W1]
    s1_D22 = s1.D[Y1,U1]

    s2_B1 = s2.B[:,W2]
    s2_B2 = s2.B[:,U2]
    s2_C1 = s2.C[Z2,:]
    s2_C2 = s2.C[Y2,:]
    s2_D11 = s2.D[Z2,W2]
    s2_D12 = s2.D[Z2,U2]
    s2_D21 = s2.D[Y2,W2]
    s2_D22 = s2.D[Y2,U2]


    if iszero(s1_D22) || iszero(s2_D22)
        A = [s1.A + α*s1_B2*s2_D22*s1_C2        α*s1_B2*s2_C2;
                 s2_B2*s1_C2            s2.A + α*s2_B2*s1_D22*s2_C2]

        B = [s1_B1 + α*s1_B2*s2_D22*s1_D21        α*s1_B2*s2_D21;
                      s2_B2*s1_D21            s2_B1 + α*s2_B2*s1_D22*s2_D21]
        C = [s1_C1 + α*s1_D12*s2_D22*s1_C2        α*s1_D12*s2_C2;
                      s2_D12*s1_C2           s2_C1 + α*s2_D12*s1_D22*s2_C2]
        D = [s1_D11 + α*s1_D12*s2_D22*s1_D21        α*s1_D12*s2_D21;
                      s2_D12*s1_D21           s2_D11 + α*s2_D12*s1_D22*s2_D21]
    else
        R1 = inv(I - α*s2_D22*s1_D22) # inv seems to be better than lu
        R2 = inv(I - α*s1_D22*s2_D22)

        A = [s1.A + α*s1_B2*R1*s2_D22*s1_C2        α*s1_B2*R1*s2_C2;
                 s2_B2*R2*s1_C2            s2.A + α*s2_B2*R2*s1_D22*s2_C2]

        B = [s1_B1 + α*s1_B2*R1*s2_D22*s1_D21        α*s1_B2*R1*s2_D21;
                     s2_B2*R2*s1_D21            s2_B1 + α*s2_B2*R2*s1_D22*s2_D21]
        C = [s1_C1 + α*s1_D12*R1*s2_D22*s1_C2        α*s1_D12*R1*s2_C2;
                     s2_D12*R2*s1_C2           s2_C1 + α*s2_D12*R2*s1_D22*s2_C2]
        D = [s1_D11 + α*s1_D12*R1*s2_D22*s1_D21        α*s1_D12*R1*s2_D21;
                     s2_D12*R2*s1_D21           s2_D11 + α*s2_D12*R2*s1_D22*s2_D21]
    end
    # Return while reorganizing into the desired order
    return StateSpace(A, isnothing(w_reorder) ? B : B[:, w_reorder],
                         isnothing(z_reorder) ? C : C[z_reorder,:],
                         isnothing(w_reorder) && isnothing(z_reorder) ? D : D[z_reorder, w_reorder])
    #return StateSpace(A, B_tmp, C_tmp, D_tmp)
end

"""
`feedback2dof(P,R,S,T)` Return `BT/(AR+ST)` where B and A are the numerator and denomenator polynomials of `P` respectively
`feedback2dof(B,A,R,S,T)` Return `BT/(AR+ST)`
"""
function feedback2dof(P::TransferFunction,R,S,T)
    !issiso(P) && error("Feedback not implemented for MIMO systems")
    tf(conv(poly2vec(numpoly(P)[1]),T),zpconv(poly2vec(denpoly(P)[1]),R,poly2vec(numpoly(P)[1]),S))
end

feedback2dof(B,A,R,S,T) = tf(conv(B,T),zpconv(A,R,B,S))


"""
    lft(G, Δ, type=:l)

Upper linear fractional transformation between systems `G` and `Δ`.
`G` must have more inputs and outputs than `Δ` has outputs and inputs.

For details, see Chapter 9.1 in
**Zhou, K. and JC Doyle**. Essentials of robust control, Prentice hall (NJ), 1998
"""
function lft(G, Δ, type=:l) # QUESTION: Or lft_u, and lft_l instead?

    if !(G.nu > Δ.ny && G.ny > Δ.nu)
        error("Must have G.nu > Δ.ny and G.ny > Δ.nu for lower/upper lft.")
    end

    if type == :l
        feedback(G, Δ, U1=G.nu-Δ.ny+1:G.nu, Y1=G.ny-Δ.nu+1:G.ny, W1=1:G.ny-Δ.nu, Z1=1:G.nu-Δ.ny, neg_feedback=false)
    elseif type == :u
        feedback(G, Δ, U1=1:Δ.ny, Y1=1:Δ.nu, W1=Δ.nu+1:G.ny, Z1=Δ.nu+1:G.ny, neg_feedback=false)
    else
        error("Invalid type of lft ($type). Specify type=:l (:u) for lower (upper) lft.")
    end
end


"""
    starprod(sys1, sys2, dimu, dimy)

Compute the Redheffer star product.

`length(U1) = length(Y2) = dimu` and `length(Y1) = length(U2) = dimy`

For details, see Chapter 9.3 in
**Zhou, K. and JC Doyle**. Essentials of robust control, Prentice hall (NJ), 1998
"""
function starprod(G1, G2, dimy::Int, dimu::Int)
    feedback(G1, G2, U1=G1.nu-dimu+1:G1.nu, Y1=G1.ny-dimy+1:G1.ny, W1=1:G1.nu-dimu, Z1=1:G1.ny-dimy,
                     U2=1:dimy, Y2=1:dimu, W2=dimy+1:G2.nu, Z2=dimu+1:G2.ny,
                     neg_feedback=false)
end
function starprod(sys1, sys2)
    nu = size(sys2, 2)
    ny = size(sys2, 1)
    starp(sys1, sys2, nu, ny)
end
