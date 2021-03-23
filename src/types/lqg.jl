
import Base.getindex

"""
    G = LQG(sys::AbstractStateSpace, Q1, Q2, R1, R2; qQ=0, qR=0, integrator=false, M = I, N = I)

Return an LQG object that describes the closed control loop around the process `sys=ss(A,B,C,D)`
where the controller is of LQG-type. The controller is specified by weight matrices `Q1,Q2`
that penalizes state deviations and control signal variance respectively, and covariance
matrices `R1,R2` which specify state drift and measurement covariance respectively.
This constructor calls [`lqr`](@ref) and [`kalman`](@ref) and forms the closed-loop system.

If `integrator=true`, the resulting controller will have integral action.
This is achieved by adding a model of a constant disturbance on the inputs to the system
described by `sys`.

`qQ` and `qR` can be set to incorporate loop transfer recovery, i.e.,
```julia
L = lqr(A, B, Q1+qQ*C'C, Q2)
K = kalman(A, C, R1+qR*B*B', R2)
```

`M` is a matrix that defines the controlled variables `z`, i.e., the variables for which you provide reference signals. If no `M` is provided, the default is to consider all state variables of the system as controlled. The definitions of `z` and `y` are given below
```
y = C*x
z = M*x
```
`size(M, 1)` determines the size of the `Q1` matrix you need to supply.

`N` is a matrix that defines how the dynamics noise `v` enters the system, i.e. If no `N` is provided, the default is to consider all state variables being affected by independent noise components. The definition of `v` is given below
```
x′ = A*x + B*u + N*v
```
`size(N, 2)` determines the size of the `R1` matrix you need to supply.

# Fields and properties
When the LQG-object is populated by the lqg-function, the following fields have been made available
- `L` is the feedback matrix, such that `A-BL` is stable. Note that the length of the state vector (and the width of L) is increased by the number of inputs if the option `integrator=true`.
- `K` is the kalman gain such that `A-KC` is stable

Several other properties of the object are accessible as properties. The available properties are
(some have many alternative names, separated with / )

- `G.cl / G.closedloop` is the closed-loop system, including observer, from reference to output, precompensated to have static gain 1 (`u = −Lx + lᵣr`).
- `G.S / G.Sin` Input sensitivity function
- `G.T / G.Tin` Input complementary sensitivity function
- `G.Sout` Output sensitivity function
- `G.Tout` Output complementary sensitivity function
- `G.CS` The transfer function from measurement noise to control signal
- `G.DS` The transfer function from input load disturbance to output
- `G.lt / G.looptransfer / G.loopgain  =  PC`
- `G.rd / G.returndifference  =  I + PC`
- `G.sr / G.stabilityrobustness  =  I + inv(PC)`
- `G.Fy / G.controller` Returns the controller as a StateSpace-system `u = L*inv(sI - A + BL + KC)*K * y`. This controller is acting on the measured signal, not the reference. The controller acting on the reference is `G.Fr`
- `G.Fr` Returns the controller from reference as a StateSpace-system. `I - L*inv(sI - A + BL + KC)*B`

It is also possible to access all fileds using the `G.symbol` syntax, the fields are `P,Q1,Q2,R1,R2,qQ,qR,sysc,L,K,integrator`

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
    L::AbstractMatrix
    K::AbstractMatrix
    M::AbstractMatrix
    N::AbstractMatrix
    syse::StateSpace
    Lr
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
    kwargs...,
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
    kwargs...,
)
    Q1 = diagm(0 => Q1)
    Q2 = diagm(0 => Q2)
    R1 = diagm(0 => R1)
    R2 = diagm(0 => R2)
    integrator ? _LQGi(sys, Q1, Q2, R1, R2, qQ, qR; kwargs...) :
    _LQG(sys, Q1, Q2, R1, R2, qQ, qR; kwargs...)
end # (2) Dispatches to final

# This function does the actual initialization in the standard case withput integrator
function _LQG(sys::LTISystem, Q1, Q2, R1, R2, qQ, qR; M = I(nstates(sys)), N = I(nstates(sys)))
    A, B, C, D = ssdata(sys)
    n = size(A, 1)
    m = size(B, 2)
    p = size(C, 1)
    size(Q1, 1) == size(M,1) || throw(ArgumentError("The size of Q1 is determined by M, not by the state."))
    size(R2, 1) == size(C,1) || throw(ArgumentError("The size of R2 is determined by C, not by the state."))
    size(R1, 1) == size(N,2) || throw(ArgumentError("The size of R1 is determined by N, not by the state."))
    L = lqr(A, B, M'Q1*M + qQ * C'C, Q2)
    #               Q                 R
    K = kalman(A, C, N*R1*N' + qR * B * B', R2)
    #                  Q                       R

    Lr = pinv(M * ((B * L - A) \ B))
    return LQG(sys, Q1, Q2, R1, R2, qQ, qR, L, K, M, N, sys, Lr)
end


# This function does the actual initialization in the integrator case
function _LQGi(sys::LTISystem, Q1, Q2, R1, R2, qQ, qR; M = I(nstates(sys)), N = nothing, ϵ=1e-3, measurement=false)
    A, B, C, D = ssdata(sys)
    n = size(A, 1)
    m = size(B, 2)
    p = size(C, 1)
    pm = size(M, 1)

    Me = [M zeros(pm, m)] # the extension is done in getproperty

    syse = add_low_frequency_disturbance(sys; ϵ, measurement)
    Ae, Be, Ce, De = ssdata(syse)
    
    size(M, 2) == n || throw(ArgumentError("The size of M does not match the size of the A-matrix, the system has $(n) states."))
    size(Q1, 1) == size(M,1) || throw(ArgumentError("The size of Q1 is determined by M, not by the state. With the current M, you need a Q1 matrix of size $(size(M,1))"))

    if N === nothing
        N = zeros(n + m, m)
        N[end-m+1:end,:] .= I(m)
        @info "Choosing an N matrix automatically" N
    end

    size(N, 1) == size(Ae, 1) || throw(ArgumentError("The size of N does not match the size of the extended A-matrix, the extended system has $(size(Ae,1)) states."))
    size(R1, 1) == size(N,2) || throw(ArgumentError("The size of R1 is determined by N, not by the state. With the current N, you need a R1 matrix of size $(size(N,2))"))
    
    T = eltype(A)    

    L = lqr(A, B, M'Q1*M + qQ * C'C, Q2)
    Le = [L I]
    K = kalman(Ae, Ce, N*R1*N' + qR * Be * Be', R2)
    Lr = pinv(M * ((B * L - A) \ B))

    LQG(sys, Q1, Q2, R1, R2, qQ, qR, Le, K, M, N, syse, Lr)
end

@deprecate getindex(G::LQG, s::Symbol) getfield(G, s)

function Base.getproperty(G::LQG, s::Symbol)
    if s ∈ (:L, :K, :Q1, :Q2, :R1, :R2, :qQ, :qR, :integrator, :P, :M, :N, :syse, :Lr)
        return getfield(G, s)
    end
    s === :A && return G.P.A
    s === :B && return G.P.B
    s === :C && return G.P.C
    s === :D && return G.P.D
    s ∈ (:sys, :P) && return getfield(G, :P)

    A = G.P.A
    B = G.P.B
    C = G.P.C
    D = G.P.D
    M = G.M
    N = G.N

    L = G.L
    K = G.K
    P = G.P
    Lr = G.Lr

    sysc = Fy = let
        Ae,Be,Ce,De = ssdata(G.syse)
        Ac = Ae - Be*L - K*Ce + K*De*L # 8.26b
        Bc = K
        Cc = L
        Dc = 0
        ss(Ac, Bc, Cc, Dc)
    end

    n = size(A, 1)
    m = size(B, 2)
    p = size(C, 1)
    pm = size(M, 1)
    pn = size(N, 2)

    # Extract interesting values
    if G.syse != G.sys 
        A, B, C, D = ssdata(G.syse)
        Me = [M zeros(pm, m)]
    else
        Me = M
    end
    PC = P * sysc # Loop gain
    if s ∈ (:cl, :closedloop, :ry) # Closed-loop system
        # Compensate for static gain, pp. 264 G.L.
        Lr = pinv(M * ((P.B * L[:, 1:n] - P.A) \ P.B))
        if any(!isfinite, Lr) || all(iszero, Lr)
            @warn "Could not compensate for static gain automatically." Lr
            Lr = 1
        end
        Acl = [A-B*L B*L; zero(A) A-K*C] # 8.28
        BLr = B * Lr
        Bcl = [BLr; zero(BLr)]
        Ccl = [Me zero(Me)]
        syscl = ss(Acl, Bcl, Ccl, 0)
        return syscl
    elseif s ∈ (:Sout, :S) # Sensitivity function
        return output_sensitivity(P, sysc)
    elseif s ∈ (:Tout, :T) # Complementary sensitivity function
        return output_comp_sensitivity(P,sysc)
    elseif s === :Sin # Sensitivity function, input
        # return feedback(ss(Matrix{numeric_type(sysc)}(I, m, m)), sysc * P)
        return input_sensitivity(P, sysc)
    elseif s === :Tin # Complementary sensitivity function, output
        return input_comp_sensitivity(P,sysc)
    elseif s === :PS # Load disturbance to output
        return P * G.Sin
    elseif s === :CS # Noise to control signal
        return sysc * G.Sout
    elseif s ∈ (:lt, :looptransfer, :loopgain)
        return PC
    elseif s ∈ (:rd, :returndifference)
        return ss(Matrix{numeric_type(PC)}(I, p, p)) + PC
    elseif s ∈ (:sr, :stabilityrobustness)
        return ss(Matrix{numeric_type(PC)}(I, p, p)) + inv(PC)
    elseif s ∈ (:Fy, :sysc, :controller)
        return Fy
    elseif s === :Fr
        Ae,Be,Ce,De = ssdata(G.syse)
        Ac = Ae - Be*L - K*Ce + K*De*L # 8.26b
        Bc = Be*G.Lr
        Cc = L
        Dc = 0
        return 1 - ss(Ac, Bc, Cc, Dc)
    end
    error("The symbol $s does not have a function associated with it.")
end

Base.propertynames(G::LQG, private::Bool = false) = (fieldnames(typeof(G))..., :Fy, :Fr, :Sin, :Sout, :Tin, :Tout, :PS, :CS, :loopgain, :returndifference, :stabilityrobustness, :cl)

Base.:(==)(G1::LQG, G2::LQG) =
    G1.K == G2.K && G1.L == G2.L && G1.P == G2.P && G1.sysc == G2.sysc



function input_sensitivity(P,C)
    T = feedback(C * P)
    ss(I(noutputs(T))) - T
end

function input_comp_sensitivity(P,C)
    T = feedback(C * P)
end

function output_sensitivity(P,C)
    PC = P*C
    S = feedback(ss(Matrix{numeric_type(PC)}(I, ninputs(PC), ninputs(PC))), PC)
    S.C .*= -1
    S.B .*= -1
    S
end

function output_comp_sensitivity(P,C)
    S = output_sensitivity(P,C)
    ss(I(noutputs(S))) - S
end


plot(G::LQG) = gangoffourplot(G)

function gangoffour(G::LQG)
    G.S, G.PS, G.CS, G.T
end

function gangoffourplot(G::LQG, args...; kwargs...)
    S,D,N,T = gangoffour(G)
    f1 = sigmaplot(S, args...; show=false, title="\$S = 1/(1+PC)\$", kwargs...)
    f2 = sigmaplot(D, args...; show=false, title="\$D = P/(1+PC)\$", kwargs...)
    f3 = sigmaplot(N, args...; show=false, title="\$N = C/(1+PC)\$", kwargs...)
    f4 = sigmaplot(T, args...; show=false, title="\$T = PC/(1+PC\$)", kwargs...)
    Plots.plot(f1,f2,f3,f4)
end



# function gangoffourplot(G::LQG, args...)
#     S,D,N,T = gangoffour(G)
#     fig = subplot(n=4,nc=2)
#     Plots.plot!(fig[1,1],sigmaplot(S, args...), title="\$S = 1/(1+PC)\$")
#     Plots.plot!(fig[1,2],sigmaplot(D, args...), title="\$D = P/(1+PC)\$")
#     Plots.plot!(fig[2,1],sigmaplot(N, args...), title="\$N = C/(1+PC)\$")
#     Plots.plot!(fig[2,2],sigmaplot(T, args...), title="\$T = PC/(1+PC\$)")
#     return fig
# end
