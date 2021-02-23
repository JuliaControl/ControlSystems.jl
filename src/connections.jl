# Model interconnections

"""
    series(sys1::LTISystem, sys2::LTISystem)

Connect systems in series, equivalent to `sys2*sys1`
"""
series(sys1::LTISystem, sys2::LTISystem) = sys2*sys1

"""
    parallel(sys1::LTISystem, sys2::LTISystem)

Connect systems in parallel, equivalent to `sys2+sys1`
"""
parallel(sys1::LTISystem, sys2::LTISystem) = sys1 + sys2

append() = LTISystem[]
"""
    append(systems::StateSpace...), append(systems::TransferFunction...)

Append systems in block diagonal form
"""
function append(systems::(ST where ST<:AbstractStateSpace)...)
    ST = promote_type(typeof.(systems)...)
    timeevol = common_timeevol(systems...)
    A = blockdiag([s.A for s in systems]...)
    B = blockdiag([s.B for s in systems]...)
    C = blockdiag([s.C for s in systems]...)
    D = blockdiag([s.D for s in systems]...)
    return ST(A, B, C, D, timeevol)
end

function append(systems::TransferFunction...)
    timeevol = common_timeevol(systems...)
    mat = blockdiag([s.matrix for s in systems]...)
    return TransferFunction(mat, timeevol)
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
    timeevol = common_timeevol(systems...)
    return ST(A, B, C, D, timeevol)
end

function Base.vcat(systems::TransferFunction...)
    # Perform checks
    nu = systems[1].nu
    if !all(s.nu == nu for s in systems)
        error("All systems must have same input dimension")
    end
    timeevol = common_timeevol(systems...)
    mat = vcat([s.matrix for s in systems]...)

    return TransferFunction(mat, timeevol)
end

Base.vcat(systems::LTISystem...) = vcat(promote(systems...)...)

function Base.hcat(systems::ST...) where ST <: AbstractStateSpace
    # Perform checks
    ny = systems[1].ny
    if !all(s.ny == ny for s in systems)
        error("All systems must have same output dimension")
    end
    timeevol = common_timeevol(systems...)
    A = blockdiag([s.A for s in systems]...)
    B = blockdiag([s.B for s in systems]...)
    C = hcat([s.C for s in systems]...)
    D = hcat([s.D for s in systems]...)

    return ST(A, B, C, D, timeevol)
end

function Base.hcat(systems::TransferFunction...)
    # Perform checks
    ny = systems[1].ny
    if !all(s.ny == ny for s in systems)
        error("All systems must have same output dimension")
    end
    timeevol = common_timeevol(systems...)
    mat = hcat([s.matrix for s in systems]...)

    return TransferFunction(mat, timeevol)
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



"""
    feedback(L)
    feedback(P1,P2)

Returns `L/(1+L)` or `P1/(1+P1*P2)`
"""
feedback(L::TransferFunction) = L/(1+L)
feedback(P1::TransferFunction, P2::TransferFunction) = P1/(1+P1*P2)

#Efficient implementations
function feedback(L::TransferFunction{<:TimeEvolution,T}) where T<:SisoRational
    if size(L) != (1,1)
        error("MIMO TransferFunction feedback isn't implemented, use L/(1+L)")
    end
    P = numpoly(L)
    Q = denpoly(L)
    tf(P, P+Q, L.timeevol)
end

function feedback(L::TransferFunction{TE, T}) where {TE<:TimeEvolution, T<:SisoZpk}
    if size(L) != (1,1)
        error("MIMO TransferFunction feedback isn't implemented, use L/(1+L)")
    end
    #Extract polynomials and create P/(P+Q)
    k = L.matrix[1].k
    denpol = numpoly(L)[1]+denpoly(L)[1]
    kden = denpol[end] # Get coeff of s^n
    # Create siso system
    sisozpk = T(L.matrix[1].z, roots(denpol), k/kden)
    return TransferFunction{TE,T}(fill(sisozpk,1,1), L.timeevol)
end

"""
    feedback(sys)
    feedback(sys1,sys2)

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
    timeevol = common_timeevol(sys1,sys2)
    !(iszero(sys1.D) || iszero(sys2.D)) && error("There cannot be a direct term (D) in both sys1 and sys2")
    A = [sys1.A+sys1.B*(-sys2.D)*sys1.C sys1.B*(-sys2.C);
         sys2.B*sys1.C  sys2.A+sys2.B*sys1.D*(-sys2.C)]
    B = [sys1.B; sys2.B*sys1.D]
    C = [sys1.C  sys1.D*(-sys2.C)]
    ss(A, B, C, sys1.D, timeevol)
end


"""
    feedback(s1::AbstractStateSpace, s2::AbstractStateSpace;
             U1=:, Y1=:, U2=:, Y2=:, W1=:, Z1=:, W2=Int[], Z2=Int[],
             Wperm=:, Zperm=:, pos_feedback::Bool=false)


`U1`, `Y1`, `U2`, `Y2` contain the indices of the signals that should be connected.
`W1`, `Z1`, `W2`, `Z2` contain the signal indices of `s1` and `s2` that should be kept.

Specify  `Wperm` and `Zperm` to reorder [w1; w2] and [z1; z2] in the resulting statespace model.

Negative feedback is the default. Specify `pos_feedback=true` for positive feedback.

See Zhou etc. for similar (somewhat less symmetric) formulas.
"""
@views function feedback(sys1::AbstractStateSpace, sys2::AbstractStateSpace;
    U1=:, Y1=:, U2=:, Y2=:, W1=:, Z1=:, W2=Int[], Z2=Int[],
    Wperm=:, Zperm=:, pos_feedback::Bool=false)

    timeevol = common_timeevol(sys1,sys2)

    if !(isa(Y1, Colon) || allunique(Y1)); @warn "Connecting single output to multiple inputs Y1=$Y1"; end
    if !(isa(Y2, Colon) || allunique(Y2)); @warn "Connecting single output to multiple inputs Y2=$Y2"; end
    if !(isa(U1, Colon) || allunique(U1)); @warn "Connecting multiple outputs to a single input U1=$U1"; end
    if !(isa(U2, Colon) || allunique(U2)); @warn "Connecting a single output to multiple inputs U2=$U2"; end

    if (U1 isa Colon ? size(sys1, 2) : length(U1)) != (Y2 isa Colon ? size(sys2, 1) : length(Y2))
        error("Lengths of U1 ($U1) and Y2 ($Y2) must be equal")
    end
    if (U2 isa Colon ? size(sys2, 2) : length(U2)) != (Y1 isa Colon ? size(sys1, 1) : length(Y1))
        error("Lengths of U1 ($U2) and Y2 ($Y1) must be equal")
    end

    α = pos_feedback ? 1 : -1 # The sign of feedback

    s1_B1 = sys1.B[:,W1]
    s1_B2 = sys1.B[:,U1]
    s1_C1 = sys1.C[Z1,:]
    s1_C2 = sys1.C[Y1,:]
    s1_D11 = sys1.D[Z1,W1]
    s1_D12 = sys1.D[Z1,U1]
    s1_D21 = sys1.D[Y1,W1]
    s1_D22 = sys1.D[Y1,U1]

    s2_B1 = sys2.B[:,W2]
    s2_B2 = sys2.B[:,U2]
    s2_C1 = sys2.C[Z2,:]
    s2_C2 = sys2.C[Y2,:]
    s2_D11 = sys2.D[Z2,W2]
    s2_D12 = sys2.D[Z2,U2]
    s2_D21 = sys2.D[Y2,W2]
    s2_D22 = sys2.D[Y2,U2]

    if iszero(s1_D22) || iszero(s2_D22)
        A = [sys1.A + α*s1_B2*s2_D22*s1_C2        α*s1_B2*s2_C2;
                 s2_B2*s1_C2            sys2.A + α*s2_B2*s1_D22*s2_C2]

        B = [s1_B1 + α*s1_B2*s2_D22*s1_D21        α*s1_B2*s2_D21;
                      s2_B2*s1_D21            s2_B1 + α*s2_B2*s1_D22*s2_D21]
        C = [s1_C1 + α*s1_D12*s2_D22*s1_C2        α*s1_D12*s2_C2;
                      s2_D12*s1_C2           s2_C1 + α*s2_D12*s1_D22*s2_C2]
        D = [s1_D11 + α*s1_D12*s2_D22*s1_D21        α*s1_D12*s2_D21;
                      s2_D12*s1_D21           s2_D11 + α*s2_D12*s1_D22*s2_D21]
    else
        # inv seems to be better than lu
        R1 = try
            inv(α*I - s2_D22*s1_D22) # slightly faster than α*inv(I - α*s2_D22*s1_D22)
        catch
            error("Ill-posed feedback interconnection,  I - α*s2_D22*s1_D22 or I - α*s2_D22*s1_D22 not invertible")
        end

        R2 = try
            inv(I - α*s1_D22*s2_D22)
        catch
            error("Ill-posed feedback interconnection,  I - α*s2_D22*s1_D22 or I - α*s2_D22*s1_D22 not invertible")
        end

        A = [sys1.A + s1_B2*R1*s2_D22*s1_C2        s1_B2*R1*s2_C2;
                 s2_B2*R2*s1_C2            sys2.A + α*s2_B2*R2*s1_D22*s2_C2]

        B = [s1_B1 + s1_B2*R1*s2_D22*s1_D21        s1_B2*R1*s2_D21;
                     s2_B2*R2*s1_D21            s2_B1 + α*s2_B2*R2*s1_D22*s2_D21]
        C = [s1_C1 + s1_D12*R1*s2_D22*s1_C2        s1_D12*R1*s2_C2;
                     s2_D12*R2*s1_C2           s2_C1 + α*s2_D12*R2*s1_D22*s2_C2]
        D = [s1_D11 + s1_D12*R1*s2_D22*s1_D21        s1_D12*R1*s2_D21;
                     s2_D12*R2*s1_D21           s2_D11 + α*s2_D12*R2*s1_D22*s2_D21]
    end

    return StateSpace(A, B[:, Wperm], C[Zperm,:], D[Zperm, Wperm], timeevol)
end


"""
    feedback2dof(P,R,S,T)
    feedback2dof(B,A,R,S,T)

- Return `BT/(AR+ST)` where B and A are the numerator and denomenator polynomials of `P` respectively
- Return `BT/(AR+ST)`
"""
function feedback2dof(P::TransferFunction,R,S,T)
    !issiso(P) && error("Feedback not implemented for MIMO systems")
    tf(conv(poly2vec(numpoly(P)[1]),T),zpconv(poly2vec(denpoly(P)[1]),R,poly2vec(numpoly(P)[1]),S), P.timeevol)
end

feedback2dof(B,A,R,S,T) = tf(conv(B,T),zpconv(A,R,B,S))

"""
    feedback2dof(P::TransferFunction, C::TransferFunction, F::TransferFunction)

Return the transfer function
`P(F+C)/(1+PC)`
which is the closed-loop system with process `P`, controller `C` and feedforward filter `F` from reference to control signal (by-passing `C`).

         +-------+
         |       |
   +----->   F   +----+
   |     |       |    |
   |     +-------+    |
   |     +-------+    |    +-------+
r  |  -  |       |    |    |       |    y
+--+----->   C   +----+---->   P   +---+-->
      |  |       |         |       |   |
      |  +-------+         +-------+   |
      |                                |
      +--------------------------------+
"""
function feedback2dof(P::TransferFunction{TE}, C::TransferFunction{TE}, F::TransferFunction{TE}) where TE
    !issiso(P) && error("Feedback not implemented for MIMO systems")
    timeevol = common_timeevol(P, C, F)
    
    Pn,Pd = numpoly(P)[], denpoly(P)[]
    Cn,Cd = numpoly(C)[], denpoly(C)[]
    Fn,Fd = numpoly(F)[], denpoly(F)[]
    den = (Cd*Pd + Pn*Cn)*Fd
    tf(Cd*Pn*Fn + Pn*Cn*Fd, den, timeevol)
end

"""
    lft(G, Δ, type=:l)

Lower and upper linear fractional transformation between systems `G` and `Δ`.

Specify `:l` lor lower LFT, and `:u` for upper LFT.

`G` must have more inputs and outputs than `Δ` has outputs and inputs.

For details, see Chapter 9.1 in
**Zhou, K. and JC Doyle**. Essentials of robust control, Prentice hall (NJ), 1998
"""
function lft(G, Δ, type=:l)

    if !(G.nu > Δ.ny && G.ny > Δ.nu)
        error("Must have G.nu > Δ.ny and G.ny > Δ.nu for lower/upper lft")
    end

    if type === :l
        feedback(G, Δ, U1=G.nu-Δ.ny+1:G.nu, Y1=G.ny-Δ.nu+1:G.ny, W1=1:G.ny-Δ.nu, Z1=1:G.nu-Δ.ny, pos_feedback=true)
    elseif type === :u
        feedback(G, Δ, U1=1:Δ.ny, Y1=1:Δ.nu, W1=Δ.nu+1:G.ny, Z1=Δ.nu+1:G.ny, pos_feedback=true)
    else
        error("Invalid type of lft ($type), specify type=:l (:u) for lower (upper) lft")
    end
end



"""
    starprod(sys1, sys2, dimu, dimy)

Compute the Redheffer star product.

`length(U1) = length(Y2) = dimu` and `length(Y1) = length(U2) = dimy`

For details, see Chapter 9.3 in
**Zhou, K. and JC Doyle**. Essentials of robust control, Prentice hall (NJ), 1998
"""
starprod(G1, G2, dimy::Int, dimu::Int) = feedback(G1, G2,
         U1=G1.nu-dimu+1:G1.nu, Y1=G1.ny-dimy+1:G1.ny, W1=1:G1.nu-dimu, Z1=1:G1.ny-dimy,
         U2=1:dimy, Y2=1:dimu, W2=dimy+1:G2.nu, Z2=dimu+1:G2.ny,
         pos_feedback=true)
starprod(sys1, sys2) = lft(sys1, sys2, :l)
