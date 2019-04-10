
import Base.getindex

"""
    G = LQG(A,B,C,D, Q1, Q2, R1, R2; qQ=0, qR=0, integrator=false)
    G = LQG(sys, args...; kwargs...)

Return an LQG object that describes the closed control loop around the process `sys=ss(A,B,C,D)`
where the controller is of LQG-type. The controller is specified by weight matrices `Q1,Q2`
that penalizes state deviations and control signal variance respectively, and covariance
matrices `R1,R2` which specify state drift and measurement covariance respectively.
This constructor calls [`lqr`](@ref) and [`kalman`](@ref) and forms the closed-loop system.

If `integrator=true`, the resulting controller will have intregral action.
This is achieved by adding a model of a constant disturbance on the inputs to the system
described by `A,B,C,D`.

`qQ` and `qR` can be set to incorporate loop transfer recovery, i.e.,
```julia
L = lqr(A, B, Q1+qQ*C'C, Q2)
K = kalman(A, C, R1+qR*B*B', R2)
```

# Fields
When the LQG-object is populated by the lqg-function, the following fields have been made available
- `L` is the feedback matrix, such that `A-BL` is stable. Note that the length of the state vector (and the width of L) is increased by the number of inputs if the option `integrator=true`.
- `K` is the kalman gain such that `A-KC` is stable
- `sysc` is a dynamical system describing the controller `u=L*inv(A-BL-KC+KDL)Ky`

# Functions
Several other properties of the object are accessible with the indexing function `getindex()`
and are called with the syntax `G[:function]`. The available functions are
(some have many alternative names, separated with / )

-`G[:cl] / G[:closedloop]` is the closed-loop system, including observer, from reference to output, precompensated to have static gain 1 (`u = −Lx + lᵣr`).
-`G[:S] / G[:Sin]` Input sensitivity function
-`G[:T] / G[:Tin]` Input complementary sensitivity function
-`G[:Sout]` Output sensitivity function
-`G[:Tout]` Output complementary sensitivity function
-`G[:CS]` The transfer function from measurement noise to control signal
-`G[:DS]` The transfer function from input load disturbance to output
-`G[:lt] / G[:looptransfer] / G[:loopgain]  =  PC`
-`G[:rd] / G[:returndifference]  =  I + PC`
-`G[:sr] / G[:stabilityrobustness]  =  I + inv(PC)`
-`G[:sysc] / G[:controller]` Returns the controller as a StateSpace-system

It is also possible to access all fileds using the `G[:symbol]` syntax, the fields are `P
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

Gcl = G[:cl]
T = G[:T]
S = G[:S]
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
    integrator::Bool
end

# Provide some constructors
function LQG(A,B,C,D,Q1::AbstractMatrix,Q2::AbstractMatrix,R1::AbstractMatrix,R2::AbstractMatrix; qQ=0, qR=0, integrator=false)
    integrator ? _LQGi(A,B,C,D,Q1,Q2,R1,R2, qQ, qR) : _LQG(A,B,C,D,Q1,Q2,R1,R2, qQ, qR)
end # (1) Dispatches to final

function LQG(A,B,C,D,Q1::AbstractVector,Q2::AbstractVector,R1::AbstractVector,R2::AbstractVector; qQ=0, qR=0, integrator=false)
    Q1 = diagm(0 => Q1)
    Q2 = diagm(0 => Q2)
    R1 = diagm(0 => R1)
    R2 = diagm(0 => R2)
    integrator ? _LQGi(A,B,C,D,Q1,Q2,R1,R2, qQ, qR) : _LQG(A,B,C,D,Q1,Q2,R1,R2, qQ, qR)
end # (2) Dispatches to final

# (3) For conveniece of sending a sys, dispatches to (1/2)
LQG(sys::LTISystem, args...; kwargs...) = LQG(sys.A,sys.B,sys.C,sys.D,args...; kwargs...)



# This function does the actual initialization in the standard case withput integrator
function _LQG(A,B,C,D, Q1, Q2, R1, R2, qQ, qR)
    n = size(A,1)
    m = size(B,2)
    p = size(C,1)
    L = lqr(A, B, Q1+qQ*C'C, Q2)
    K = kalman(A, C, R1+qR*B*B', R2)

    # Controller system
    Ac=A-B*L-K*C+K*D*L
    Bc=K
    Cc=L
    Dc=zero(D')
    sysc = ss(Ac,Bc,Cc,Dc)

    return LQG(ss(A,B,C,D),Q1,Q2,R1,R2, qQ, qR, sysc, L, K, false)
end


# This function does the actual initialization in the integrator case
function _LQGi(A,B,C,D, Q1, Q2, R1, R2, qQ, qR)
    n = size(A,1)
    m = size(B,2)
    p = size(C,1)

    # Augment with disturbance model
    Ae = [A B; zeros(m,n+m)]
    Be = [B;zeros(m,m)]
    Ce = [C zeros(p,m)]
    De = D

    L = lqr(A, B, Q1+qQ*C'C, Q2)
    Le = [L I]
    K = kalman(Ae, Ce, R1+qR*Be*Be', R2)

    # Controller system
    Ac=Ae-Be*Le-K*Ce+K*De*Le
    Bc=K
    Cc=Le
    Dc=zero(D')
    sysc = ss(Ac,Bc,Cc,Dc)

    LQG(ss(A,B,C,D),Q1,Q2,R1,R2, qQ, qR, sysc, Le, K, true)
end


function Base.getindex(G::LQG, s)
    s == :A && return G.P.A
    s == :B && return G.P.B
    s == :C && return G.P.C
    s == :D && return G.P.D
    s == :L  && return G.L
    s == :K  && return G.K
    s == :Q1 && return G.Q1
    s == :Q2 && return G.Q2
    s == :R1 && return G.R1
    s == :R2 && return G.R2
    s == :qQ && return G.qQ
    s == :qR && return G.qR
    s ∈ [:sys, :P]  && return G.P
    s ∈ [:sysc, :controller] && return G.sysc
    s == :integrator && return G.integrator

    A = G.P.A
    B = G.P.B
    C = G.P.C
    D = G.P.D

    L = G.L
    K = G.K
    P = G.P
    sysc = G.sysc

    n = size(A,1)
    m = size(B,2)
    p = size(C,1)

    # Extract interesting values
    if G.integrator # Augment with disturbance model
        A = [A B; zeros(m,n+m)]
        B = [B;zeros(m,m)]
        C = [C zeros(p,m)]
        D = D
    end

    PC = P*sysc # Loop gain

    if s ∈ [:cl, :closedloop, :ry] # Closed-loop system
        Acl = [A-B*L B*L; zero(A) A-K*C]
        Bcl = [B; zero(B)]
        Ccl = [C zero(C)]
        Bcl = Bcl/(P.C*inv(P.B*L[:,1:n]-P.A)*P.B) # B*lᵣ # Always normalized with nominal plant static gain
        return syscl = ss(Acl,Bcl,Ccl,0)
    elseif s ∈ [:Sin, :S] # Sensitivity function
        return feedback(ss(Matrix{numeric_type(PC)}(I, m, m)),PC)
    elseif s ∈ [:Tin, :T] # Complementary sensitivity function
        return feedback(PC)
    elseif s == :Sout # Sensitivity function, output
        return feedback(ss(Matrix{numeric_type(sys_c)}(I, m, m)),sysc*P)
    elseif s == :Tout # Complementary sensitivity function, output
        return feedback(sysc*P)
    elseif s == :PS # Load disturbance to output
        return P*G[:S]
    elseif s == :CS # Noise to control signal
        return sysc*G[:S]
    elseif s ∈ [:lt, :looptransfer, :loopgain]
        return PC
    elseif s ∈ [:rd, :returndifference]
        return  ss(Matrix{numeric_type(PC)}(I, p, p)) + PC
    elseif s ∈ [:sr, :stabilityrobustness]
        return  ss(Matrix{numeric_type(PC)}(I, p, p)) + inv(PC)
    end
    error("The symbol $s does not have a function associated with it.")
end

Base.:(==)(G1::LQG, G2::LQG) = G1.K == G2.K && G1.L == G2.L && G1.P == G2.P && G1.sysc == G2.sysc


plot(G::LQG) = gangoffourplot(G)
function gangoffour(G::LQG)
    return G[:S], G[:PS], G[:CS], G[:T]
end

function gangoffourplot(G::LQG; kwargs...)
    S,D,N,T = gangoffour(G)
    f1 = sigmaplot(S, show=false, kwargs...); Plots.plot!(title="\$S = 1/(1+PC)\$")
    f2 = sigmaplot(D, show=false, kwargs...); Plots.plot!(title="\$D = P/(1+PC)\$")
    f3 = sigmaplot(N, show=false, kwargs...); Plots.plot!(title="\$N = C/(1+PC)\$")
    f4 = sigmaplot(T, show=false, kwargs...); Plots.plot!(title="\$T = PC/(1+PC\$)")
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
