
"""
    G = LQG(sys::AbstractStateSpace, Q1, Q2, R1, R2; qQ=0, qR=0, integrator=false, M = sys.C)

Return an LQG object that describes the closed control loop around the process `sys=ss(A,B,C,D)`
where the controller is of LQG-type. The controller is specified by weight matrices `Q1,Q2`
that penalizes state deviations and control signal variance respectively, and covariance
matrices `R1,R2` which specify state drift and measurement covariance respectively.
This constructor calls [`lqr`](@ref) and [`kalman`](@ref) and forms the closed-loop system.

If `integrator=true`, the resulting controller will have intregral action.
This is achieved by adding a model of a constant disturbance on the inputs to the system
described by `sys`.

`qQ` and `qR` can be set to incorporate loop transfer recovery, i.e.,
```julia
L = lqr(A, B, Q1+qQ*C'C, Q2)
K = kalman(A, C, R1+qR*B*B', R2)
```

`M` is a matrix that defines the controlled variables `z`, if none is provided, the default is to consider all measured outputs `y` of the system as controlled. The definitions of `z` and `y` are given below
```
y = C*x
z = M*x
```

# Fields
When the LQG-object is populated by the lqg-function, the following fields have been made available
- `L` is the feedback matrix, such that `A-BL` is stable. Note that the length of the state vector (and the width of L) is increased by the number of inputs if the option `integrator=true`.
- `K` is the kalman gain such that `A-KC` is stable
- `sysc` is a dynamical system describing the controller `u=L*inv(A-BL-KC+KDL)Ky`

Several other properties of the object are accessible as properties. The available properties are
(some have many alternative names, separated with / )

-`G.cl / G.closedloop` is the closed-loop system, including observer, from reference to output, precompensated to have static gain 1 (`u = −Lx + lᵣr`).
-`G.S / G.Sin` Input sensitivity function
-`G.T / G.Tin` Input complementary sensitivity function
-`G.Sout` Output sensitivity function
-`G.Tout` Output complementary sensitivity function
-`G.CS` The transfer function from measurement noise to control signal
-`G.DS` The transfer function from input load disturbance to output
-`G.lt / G.looptransfer / G.loopgain  =  PC`
-`G.rd / G.returndifference  =  I + PC`
-`G.sr / G.stabilityrobustness  =  I + inv(PC)`
-`G.sysc / G.controller` Returns the controller as a StateSpace-system

It is also possible to access all fileds using the `G.symbol` syntax, the fields are `P
,Q1,Q2,R1,R2,qQ,qR,sysc,L,K,integrator`

# Example

```julia
s = tf("s")
P = [1/(s+1) 2/(s+2); 1/(s+3) 1/(s-1)]
sys = ss(P)
eye(n) = Matrix{Float64}(I,n,n) # For convinience

qQ = 1
qR = 1
Q1 = 10eye(4)
Q2 = 1eye(2)
R1 = 1eye(6)
R2 = 1eye(2)

G = LQG(sys, Q1, Q2, R1, R2, qQ=qQ, qR=qR, integrator=true)

Gcl = G.cl
T = G.T
S = G.S
sigmaplot([S,T],exp10.(range(-3, stop=3, length=1000)))
stepplot(Gcl)
```

"""
struct LQG
    P::StateSpace
    Q1::AbstractMatrix
    Q2::AbstractMatrix
    R1::AbstractMatrix
    R2::AbstractMatrix
    qQ::Real
    qR::Real
    sysc::LTISystem
    L::AbstractMatrix
    K::AbstractMatrix
    M::AbstractMatrix
    integrator::Bool
end

# Provide some constructors
function LQG(
    sys::LTISystem,
    Q1::AbstractMatrix,
    Q2::AbstractMatrix,
    R1::AbstractMatrix,
    R2::AbstractMatrix;
    qQ = 0,
    qR = 0,
    integrator = false,
    kwargs...
)
    integrator ? _LQGi(sys, Q1, Q2, R1, R2, qQ, qR; kwargs...) :
    _LQG(sys, Q1, Q2, R1, R2, qQ, qR; kwargs...)
end # (1) Dispatches to final

function LQG(
    sys::LTISystem,
    Q1::AbstractVector,
    Q2::AbstractVector,
    R1::AbstractVector,
    R2::AbstractVector;
    qQ = 0,
    qR = 0,
    integrator = false,
    kwargs...
)
    Q1 = diagm(0 => Q1)
    Q2 = diagm(0 => Q2)
    R1 = diagm(0 => R1)
    R2 = diagm(0 => R2)
    integrator ? _LQGi(sys, Q1, Q2, R1, R2, qQ, qR; kwargs...) :
    _LQG(sys, Q1, Q2, R1, R2, qQ, qR; kwargs...)
end # (2) Dispatches to final

# This function does the actual initialization in the standard case withput integrator
function _LQG(sys::LTISystem, Q1, Q2, R1, R2, qQ, qR; M = sys.C)
    A, B, C, D = ssdata(sys)
    n = size(A, 1)
    m = size(B, 2)
    p = size(C, 1)
    L = lqr(A, B, Q1 + qQ * C'C, Q2)
    K = kalman(A, C, R1 + qR * B * B', R2)

    # Controller system
    Ac = A - B*L - K*C + K*D*L
    Bc = K
    Cc = L
    Dc = zero(D')
    sysc = ss(Ac, Bc, Cc, Dc)

    return LQG(ss(A, B, C, D), Q1, Q2, R1, R2, qQ, qR, sysc, L, K, M, false)
end


# This function does the actual initialization in the integrator case
function _LQGi(sys::LTISystem, Q1, Q2, R1, R2, qQ, qR; M = sys.C)
    A, B, C, D = ssdata(sys)
    n = size(A, 1)
    m = size(B, 2)
    p = size(C, 1)

    # Augment with disturbance model
    Ae = [A B; zeros(m, n + m)]
    Be = [B; zeros(m, m)]
    Ce = [C zeros(p, m)]
    De = D


    L = lqr(A, B, Q1 + qQ * C'C, Q2)
    Le = [L I]
    K = kalman(Ae, Ce, R1 + qR * Be * Be', R2)

    # Controller system
    Ac = Ae - Be*Le - K*Ce + K*De*Le
    Bc = K
    Cc = Le
    Dc = zero(D')
    sysc = ss(Ac, Bc, Cc, Dc)

    LQG(ss(A, B, C, D), Q1, Q2, R1, R2, qQ, qR, sysc, Le, K, M, true)
end

@deprecate getindex(G::LQG, s::Symbol) getfield(G, s)

function Base.getproperty(G::LQG, s::Symbol)
    if s ∈ (:L, :K, :Q1, :Q2, :R1, :R2, :qQ, :qR, :integrator, :P, :M)
        return getfield(G, s)
    end
    s === :A && return G.P.A
    s === :B && return G.P.B
    s === :C && return G.P.C
    s === :D && return G.P.D
    s ∈ (:sys, :P) && return getfield(G, :P)
    s ∈ (:sysc, :controller) && return getfield(G, :sysc)

    A = G.P.A
    B = G.P.B
    C = G.P.C
    D = G.P.D
    M = G.M

    L = G.L
    K = G.K
    P = G.P
    sysc = G.sysc

    n = size(A, 1)
    m = size(B, 2)
    p = size(C, 1)
    pm = size(M, 1)

    # Extract interesting values
    if G.integrator # Augment with disturbance model
        A = [A B; zeros(m, n + m)]
        B = [B; zeros(m, m)]
        C = [C zeros(p, m)]
        M = [M zeros(pm, m)]
        D = D
    end

    PC = P * sysc # Loop gain

    # Compensate for static gain, pp. 264 G.L.
    @show svdvals(P.B*L[:,1:n]-P.A)
    @show dcg = (P.C * inv(P.B*L[:,1:n]-P.A) * P.B)
    Acl = [A-B*L B*L; zero(A) A-K*C]
    Bcl = [B/dcg; zero(B)]
    Ccl = [M zero(M)]
    # rank(dcg) == size(A,1) && (Bcl = Bcl / dcg) # B*lᵣ # Always normalized with nominal plant static gain
    syscl = ss(Acl, Bcl, Ccl, 0)
    if s ∈ (:cl, :closedloop, :ry) # Closed-loop system
        # return ss(A-B*L, B/dcg, M, 0)
        return syscl
    elseif s ∈ (:Sin, :S) # Sensitivity function
        return feedback(ss(Matrix{numeric_type(PC)}(I, m, m)), PC)
    elseif s ∈ (:Tin, :T) # Complementary sensitivity function
        return feedback(PC)
        # return ss(Acl, I(size(Acl,1)), Ccl, 0)[1,2]
    elseif s === :Sout # Sensitivity function, output
        return feedback(ss(Matrix{numeric_type(sysc)}(I, m, m)), sysc * P)
    elseif s === :Tout # Complementary sensitivity function, output
        return feedback(sysc * P)
    elseif s === :PS # Load disturbance to output
        return P * G.S
    elseif s === :CS # Noise to control signal
        return sysc * G.S
    elseif s ∈ (:lt, :looptransfer, :loopgain)
        return PC
    elseif s ∈ (:rd, :returndifference)
        return ss(Matrix{numeric_type(PC)}(I, p, p)) + PC
    elseif s ∈ (:sr, :stabilityrobustness)
        return ss(Matrix{numeric_type(PC)}(I, p, p)) + inv(PC)
    end
    error("The symbol $s does not have a function associated with it.")
end

Base.:(==)(G1::LQG, G2::LQG) =
    G1.K == G2.K && G1.L == G2.L && G1.P == G2.P && G1.sysc == G2.sysc



plot(G::LQG) = gangoffourplot(G)
function gangoffour(G::LQG)
    return G.S, G.PS, G.CS, G.T
end

function gangoffourplot(G::LQG; kwargs...)
    S, D, N, T = gangoffour(G)
    f1 = sigmaplot(S, show = false, kwargs...)
    Plots.plot!(title = "\$S = 1/(1+PC)\$")
    f2 = sigmaplot(D, show = false, kwargs...)
    Plots.plot!(title = "\$D = P/(1+PC)\$")
    f3 = sigmaplot(N, show = false, kwargs...)
    Plots.plot!(title = "\$N = C/(1+PC)\$")
    f4 = sigmaplot(T, show = false, kwargs...)
    Plots.plot!(title = "\$T = PC/(1+PC\$)")
    Plots.plot(f1, f2, f3, f4)
end
