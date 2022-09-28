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
    A = blockdiag(s.A for s in systems)
    B = blockdiag(s.B for s in systems)
    C = blockdiag(s.C for s in systems)
    D = blockdiag(s.D for s in systems)
    return ST(A, B, C, D, timeevol)
end

function append(systems::TransferFunction...)
    timeevol = common_timeevol(systems...)
    mat = blockdiag(s.matrix for s in systems)
    return TransferFunction(mat, timeevol)
end

append(systems::LTISystem...) = append(promote(systems...))

append(systems::Union{<:Tuple, <:Base.Generator}) = append(systems...)

"""
    array2mimo(M::AbstractArray{<:LTISystem})

Take an array of `LTISystem`s and create a single MIMO system.
"""
function array2mimo(M::AbstractArray{<:LTISystem})
    rows = map(axes(M, 1)) do row
        hcat(M[row, :]...)
    end
    vcat(rows...)
end

function Base.vcat(systems::LFTT...) where LFTT <: LFTSystem
    P = vcat_1([sys.P for sys in systems]...) # See PartitionedStateSpace
    f = reduce(vcat, feedback_channel(sys) for sys in systems)
    return LFTT(P, f)
end

function Base.hcat(systems::LFTT...) where LFTT <: LFTSystem
    P = hcat_1([sys.P for sys in systems]...)  # See PartitionedStateSpace
    f = reduce(vcat, feedback_channel(sys) for sys in systems)
    return LFTT(P, f)
end


function Base.vcat(systems::ST...) where ST <: AbstractStateSpace
    # Perform checks
    nu = systems[1].nu
    if !all(s.nu == nu for s in systems)
        error("All systems must have same input dimension")
    end
    A = blockdiag(s.A for s in systems)
    B = reduce(vcat, s.B for s in systems)
    C = blockdiag(s.C for s in systems)
    D = reduce(vcat, s.D for s in systems)
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
    mat = reduce(vcat, s.matrix for s in systems)

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
    A = blockdiag(s.A for s in systems)
    B = blockdiag(s.B for s in systems)
    C = reduce(hcat, s.C for s in systems)
    D = reduce(hcat, s.D for s in systems)

    return ST(A, B, C, D, timeevol)
end

function Base.hcat(systems::TransferFunction...)
    # Perform checks
    ny = systems[1].ny
    if !all(s.ny == ny for s in systems)
        error("All systems must have same output dimension")
    end
    timeevol = common_timeevol(systems...)
    mat = reduce(hcat, s.matrix for s in systems)

    return TransferFunction(mat, timeevol)
end

Base.hcat(systems::LTISystem...) = reduce(hcat, promote(systems...))

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

"""
    add_input(sys::AbstractStateSpace, B2::AbstractArray, D2 = 0)

Add inputs to `sys` by forming
```math
x' = Ax + [B \\; B_2]u
y  = Cx + [D \\; D_2]u
```
If `B2` is an integer it will be interpreted as an index and an input matrix containing a single 1 at the specified index will be used.

Example:
The following example forms an innovation model that takes innovations as inputs
```julia
G   = ssrand(2,2,3, Ts=1)
K   = kalman(G, I(G.nx), I(G.ny))
sys = add_input(G, K)
```
"""
function add_input(sys::AbstractStateSpace, B2::AbstractArray, D2=0)
    T = promote_type(numeric_type(sys), eltype(B2), eltype(D2))
    A,B,C,D = ssdata(sys)
    D3 = D2 == 0 ? zeros(T, sys.ny, size(B2, 2)) : D2
    basetype(sys)(A, [B B2], C, [D D3], sys.timeevol)
end

function add_input(sys::AbstractStateSpace, b2::Integer)
    T = promote_type(numeric_type(sys), eltype(b2))
    B2 = zeros(T, sys.nx, 1)
    B2[b2] = 1
    add_input(sys, B2)
end

"""
    add_output(sys::AbstractStateSpace, C2::AbstractArray, D2 = 0)

Add outputs to `sys` by forming
```math
x' = Ax + Bu
y  = [C; C_2]x + [D; D_2]u
```
If `C2` is an integer it will be interpreted as an index and an output matrix containing a single 1 at the specified index will be used.
"""
function add_output(sys::AbstractStateSpace, C2::AbstractArray, D2=0)
    T = promote_type(numeric_type(sys), eltype(C2), eltype(D2))
    A,B,C,D = ssdata(sys)
    D3 = D2 == 0 ? zeros(T, size(C2, 1), sys.nu) : D2
    basetype(sys)(A, B, [C; C2], [D; D3], sys.timeevol)
end

function add_output(sys::AbstractStateSpace, c2::Integer)
    T = promote_type(numeric_type(sys), eltype(c2))
    C2 = zeros(T, 1, sys.nx)
    C2[c2] = 1
    add_output(sys, C2)
end

# Catch special cases where inv(sys) might not be possible after promotion, like improper tf
function /(sys1::Union{StateSpace,AbstractStateSpace}, sys2::LTISystem)
    sys1new, sys2new = promote(sys1, 1/sys2)
    return sys1new*sys2new
end

@static if VERSION >= v"1.8.0-beta1"
    blockdiag(anything...) = cat(anything..., dims=Val((1,2)))
    blockdiag(anything::Union{<:Tuple, <:Base.Generator}) = cat(anything..., dims=Val((1,2)))
else
    blockdiag(anything...) = cat(anything..., dims=(1,2))
    blockdiag(anything::Union{<:Tuple, <:Base.Generator}) = cat(anything..., dims=(1,2))
end


"""
    feedback(sys)
    feedback(sys1, sys2)

For a general LTI-system, `feedback` forms the negative feedback interconnection
```julia
>-+ sys1 +-->
  |      |
 (-)sys2 +
```
If no second system is given, negative identity feedback is assumed
"""
feedback(L::TransferFunction) = L/(1+L)
feedback(P1::TransferFunction, P2::TransferFunction) = P1/(1+P1*P2)

function feedback(G1::TransferFunction{<:TimeEvolution,<:SisoRational}, G2::TransferFunction{<:TimeEvolution,<:SisoRational})
    if !issiso(G1) || !issiso(G2)
        error("MIMO TransferFunction feedback isn't implemented.")
    end
    G1num = numpoly(G1)[]
    G1den = denpoly(G1)[]
    G2num = numpoly(G2)[]
    G2den = denpoly(G2)[]
    
    timeevol = common_timeevol(G1, G2)
    tf(G1num*G2den, G1num*G2num + G1den*G2den, timeevol)
end

#Efficient implementations
function feedback(L::TransferFunction{<:TimeEvolution,T}) where T<:SisoRational
    if size(L) != (1,1)
        error("MIMO TransferFunction feedback isn't implemented, use L/(1+L)")
    end
    P = numpoly(L)
    Q = denpoly(L)
    tf(P, P+Q, timeevol(L))
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
    return TransferFunction{TE,T}(fill(sisozpk,1,1), timeevol(L))
end

function feedback(sys::Union{AbstractStateSpace, LFTSystem})
    ninputs(sys) != noutputs(sys) && error("Use feedback(sys1, sys2) if number of inputs != outputs")
    feedback(sys,ss(Matrix{numeric_type(sys)}(I,size(sys)...), timeevol(sys)))
end

"""
    feedback(sys1::AbstractStateSpace, sys2::AbstractStateSpace;
             U1=:, Y1=:, U2=:, Y2=:, W1=:, Z1=:, W2=Int[], Z2=Int[],
             Wperm=:, Zperm=:, pos_feedback::Bool=false)

*Basic use* `feedback(sys1, sys2)` forms the feedback interconnection
```julia
           ┌──────────────┐
◄──────────┤     sys1     │◄──── Σ ◄──────
    │      │              │      │
    │      └──────────────┘      -1
    │                            |
    │      ┌──────────────┐      │
    └─────►│     sys2     ├──────┘
           │              │
           └──────────────┘
```
*Advanced use*
`feedback` also supports more flexible use according to the figure below
```julia
              ┌──────────────┐
      z1◄─────┤     sys1     │◄──────w1
 ┌─── y1◄─────┤              │◄──────u1 ◄─┐
 │            └──────────────┘            │
 │                                        α
 │            ┌──────────────┐            │
 └──► u2─────►│     sys2     ├───────►y2──┘
      w2─────►│              ├───────►z2
              └──────────────┘
```
`U1`, `W1` specifies the indices of the input signals of `sys1` corresponding to `u1` and `w1`
`Y1`, `Z1` specifies the indices of the output signals of `sys1` corresponding to `y1` and `z1`
`U2`, `W2`, `Y2`, `Z2` specifies the corresponding signals of `sys2` 

Specify  `Wperm` and `Zperm` to reorder the inputs (corresponding to [w1; w2])
and outputs (corresponding to [z1; z2]) in the resulting statespace model.

Negative feedback (α = -1) is the default. Specify `pos_feedback=true` for positive feedback (α = 1).

See also `lft`, `starprod`, `sensitivity`, `input_sensitivity`, `output_sensitivity`, `comp_sensitivity`, `input_comp_sensitivity`, `output_comp_sensitivity`, `G_PS`, `G_CS`.

See Zhou, Doyle, Glover (1996) for similar (somewhat less symmetric) formulas.
"""
function feedback(sys1::AbstractStateSpace, sys2::AbstractStateSpace;
    U1=:, Y1=:, U2=:, Y2=:, W1=:, Z1=:, W2=Int[], Z2=Int[],
    Wperm=:, Zperm=:, pos_feedback::Bool=false)

    timeevol = common_timeevol(sys1,sys2)
    T = Base.promote_type(numeric_type(sys1), numeric_type(sys2))

    if !(isa(Y1, Colon) || allunique(Y1)); @warn "Connecting single output to multiple inputs Y1=$Y1"; end
    if !(isa(Y2, Colon) || allunique(Y2)); @warn "Connecting single output to multiple inputs Y2=$Y2"; end
    if !(isa(U1, Colon) || allunique(U1)); @warn "Connecting multiple outputs to a single input U1=$U1"; end
    if !(isa(U2, Colon) || allunique(U2)); @warn "Connecting a single output to multiple inputs U2=$U2"; end

    if (U1 isa Colon ? size(sys1, 2) : length(U1)) != (Y2 isa Colon ? size(sys2, 1) : length(Y2))
        error("Lengths of U1 ($U1) and Y2 ($Y2) must be equal")
    end
    if (U2 isa Colon ? size(sys2, 2) : length(U2)) != (Y1 isa Colon ? size(sys1, 1) : length(Y1))
        error("Lengths of U2 ($U2) and Y1 ($Y1) must be equal")
    end

    α = pos_feedback ? 1 : -1 # The sign of feedback

    @views begin
        s1_B2 = U1 isa Colon ? sys1.B : sys1.B[:,U1]
        s1_C2 = Y1 isa Colon ? sys1.C : sys1.C[Y1,:]
        s1_D12 = Z1 isa Colon && U1 isa Colon ? sys1.D : sys1.D[Z1,U1]
        s1_D21 = Y1 isa Colon && W1 isa Colon ? sys1.D : sys1.D[Y1,W1]
        s1_D22 = Y1 isa Colon && U1 isa Colon ? sys1.D : sys1.D[Y1,U1]
        
        s2_B2 = sys2.B[:,U2]
        s2_C2 = sys2.C[Y2,:]
        s2_D22 = sys2.D[Y2,U2]
    end
    # These are deliberate copies instead of views for two reasons.
    # 1) Many of them are overwritten in-place below
    # 2) Some that are not overwritten are still copied since it greatly reduces compile time. These are the ones that are by default indexed by an empty array (W2/Z2), implying that the copy will be an empty array as well.
    s1_D11 = sys1.D[Z1,W1]
    s1_C1 = sys1.C[Z1,:]
    s1_B1 = sys1.B[:,W1]
    s2_B1 = sys2.B[:,W2]
    s2_C1 = sys2.C[Z2,:]
    s2_D11 = sys2.D[Z2,W2]
    s2_D21 = sys2.D[Y2,W2]
    s2_D12 = sys2.D[Z2,U2]

    αs2_D12 = α*s2_D12
    s2_D22s1_C2 = s2_D22*s1_C2
    αs2_B2 = α*s2_B2
    s1_D22s2_C2 = s1_D22*s2_C2
    s1_D22s2_D21 = s1_D22*s2_D21

    if iszero(s1_D22) || iszero(s2_D22)
        αs1_D12 = α*s1_D12
        A11 = mul!(mutable(copy(sys1.A), T), s1_B2, s2_D22s1_C2, α, 1)
        A22 = mul!(mutable(copy(sys2.A), T), αs2_B2, s1_D22s2_C2, 1, 1)
        C11 = mul!(mutable(s1_C1, T), αs1_D12, s2_D22s1_C2, 1, 1)
        C22 = mul!(mutable(s2_C1, T), αs2_D12, s1_D22s2_C2, 1, 1)
        B11 = mul!(mutable(s1_B1, T), s1_B2, s2_D22*s1_D21, α, 1)
        B22 = mul!(mutable(s2_B1, T), αs2_B2, s1_D22s2_D21, 1, 1)
        D22 = mul!(mutable(s2_D11, T), αs2_D12, s1_D22s2_D21, 1, 1)
        D11 = mul!(mutable(s1_D11, T), αs1_D12, s2_D22*s1_D21, 1, 1)
        A = [A11        ((s1_B2*s2_C2) .*= α);
            s2_B2*s1_C2            A22]

        B = [B11        ((s1_B2*s2_D21) .*= α);
            s2_B2*s1_D21            B22]
        C = [C11        αs1_D12*s2_C2;
                      s2_D12*s1_C2           C22]
        D = [D11        αs1_D12*s2_D21;
                      s2_D12*s1_D21           D22]
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

        s2_B2R2 = s2_B2*R2
        s1_D12R1 = s1_D12*R1
        s1_B2R1 = s1_B2*R1
        s2_D22s1_D21 = s2_D22*s1_D21
        αs2_D12R2 = αs2_D12*R2
        αs2_B2R2 = α*s2_B2R2
        A = [sys1.A + s1_B2R1*s2_D22s1_C2        s1_B2R1*s2_C2;
                 s2_B2R2*s1_C2            sys2.A + αs2_B2R2*s1_D22s2_C2]

        B = [s1_B1 + s1_B2R1*s2_D22s1_D21        s1_B2R1*s2_D21;
                     s2_B2R2*s1_D21            s2_B1 + αs2_B2R2*s1_D22s2_D21]
        C = [s1_C1 + s1_D12R1*s2_D22s1_C2        s1_D12R1*s2_C2;
                     s2_D12*R2*s1_C2           s2_C1 + αs2_D12R2*s1_D22s2_C2]
        D = [s1_D11 + s1_D12R1*s2_D22s1_D21        s1_D12R1*s2_D21;
                     s2_D12*R2*s1_D21           s2_D11 + αs2_D12R2*s1_D22s2_D21]
    end

    return StateSpace(A, B[:, Wperm], C[Zperm,:], D[Zperm, Wperm], timeevol)
end

mutable(x::AbstractArray, ::Type{T}) where T = convert(Matrix{T}, x)
mutable(x::StaticArray, ::Type{T}) where T = Matrix{T}(x)


"""
    feedback2dof(P,R,S,T)
    feedback2dof(B,A,R,S,T)

- Return `BT/(AR+ST)` where B and A are the numerator and denomenator polynomials of `P` respectively
- Return `BT/(AR+ST)`
"""
function feedback2dof(P::TransferFunction,R,S,T)
    !issiso(P) && error("Feedback not implemented for MIMO systems")
    tf(conv(poly2vec(numpoly(P)[1]),T),zpconv(poly2vec(denpoly(P)[1]),R,poly2vec(numpoly(P)[1]),S), timeevol(P))
end

feedback2dof(B,A,R,S,T) = tf(conv(B,T),zpconv(A,R,B,S))

"""
    feedback2dof(P::TransferFunction, C::TransferFunction, F::TransferFunction)

Return the transfer function
`P(F+C)/(1+PC)`
which is the closed-loop system with process `P`, controller `C` and feedforward filter `F` from reference to control signal (by-passing `C`).
```
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
```
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
        feedback(G, Δ, U1=G.nu-Δ.ny+1:G.nu, Y1=G.ny-Δ.nu+1:G.ny, W1=1:G.nu-Δ.ny, Z1=1:G.ny-Δ.nu, pos_feedback=true)
    elseif type === :u
        feedback(G, Δ, U1=1:Δ.ny, Y1=1:Δ.nu, W1=Δ.ny+1:G.nu, Z1=Δ.nu+1:G.ny, pos_feedback=true)
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
