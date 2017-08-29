
import Base.getindex

"""
Initializes an LQG object. This object can be used to contruct an LQG controller.

# Constructors

`LQG(P,Q1,Q2,R1,R2, qQ, qR, sysc, L, K, integrator)` This is mainly used inside the toolbox

`LQG(P,Q1,Q2,R1,R2; qQ=0, qR=0, integrator=false)` Supply all weighting and covariance matrices as well as optional loop-transfer recovery gains. The boolean `integrator` indicates whether or not the system will be augmented with constant input disturbances to introduce integral action.

`LQG(P, qQ, qR; integrator=false)`  if called like this, where `qQ,qR` are scalars, the weighting matrices will be
`Q1 = qQ² ⋅ C'C`
`Q2 = I`
`R1 = qR² ⋅ B B'`
`R2 = I`

# Fields

When the LQG-object is populated by the lqg-function, the following fields have been made available

`L` is the feedback matrix, such that `A-BL` is stable. Note that the length of the state vector (and the width of L) is increased by the number of inputs if the option `integrator=true`.

`K` is the kalman gain such that `A-KC` is stable

`sysc` is a dynamical system describing the controller `u=L*inv(A-BL-KC+KDL)Ky`

# Functions

Several other properties of the object are accessible with the indexing function `getindex()` and are called with the syntax `G[:function]`. The available functions are (some have many alternative names, separated with / )

`G[:cl] / G[:closedloop]` is the closed-loop system, including observer, from reference to output, precompensated to have static gain 1 (`u = −Lx + lᵣr`).

`G[:S] / G[:Sin]` Input sensitivity function

`G[:T] / G[:Tin]` Input complementary sensitivity function

`G[:Sout]` Output sensitivity function

`G[:Tout]` Output complementary sensitivity function

`G[:CS]` The transfer function from measurement noise to control signal

`G[:DS]` The transfer function from input load disturbance to output

`G[:lt] / G[:looptransfer] / G[:loopgain]  =  PC`

`G[:rd] / G[:returndifference]  =  I + PC`

`G[:sr] / G[:stabilityrobustness]  =  I + inv(PC)`

`G[:sysc] / G[:controller]` Returns the controller as a StateSpace-system

It is also possible to access all fileds using the `G[:symbol]` syntax, the fields are `P
,Q1,Q2,R1,R2,qQ,qR,sysc,L,K,integrator`

# Usage example
## Specifying everything
```julia
qQ = 1
qR = 1
Q1 = 10eye(4)
Q2 = 1eye(2)
R1 = 1eye(6)
R2 = 1eye(2)
Ginit = LQG(sys, Q1, Q2, R1, R2, qQ=qQ, qR=qR, integrator=true)
G = lqg(Ginit)

Gcl = G[:cl]
T = G[:T]
S = G[:S]
sigmaplot([S,T],logspace(-3,3,1000))
stepplot(Gcl)
```

## Lazy method
```julia
qQ = 10
qR = 10
Ginit = LQG(qQ, qR, integrator=true)
G = lqg(Ginit)

Gcl = G[:cl]
T = G[:T]
S = G[:S]
sigmaplot([S,T],logspace(-3,3,1000))
stepplot(Gcl)
```
"""
type LQG
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

LQG(P,Q1::AbstractMatrix,Q2::AbstractMatrix,R1::AbstractMatrix,R2::AbstractMatrix; qQ=0, qR=0, sysc=ss(0), L=Matrix(), K=Matrix(), integrator=false) =
LQG(P,Q1,Q2,R1,R2, qQ, qR, sysc, L, K, integrator)

function LQG(P, qQ::Real, qR::Real; sysc=ss(0), L=Matrix(), K=Matrix(), integrator=false)
    m = size(P.B,2)
    p = size(P.C,1)
    n = size(P.A,1)
    nr = n + integrator ? m : 0
    LQG(P,0eye(n),eye(m),0eye(nr),eye(p), qQ, qR, sysc, L, K, integrator)
end

function LQG(P,Q1::AbstractVector,Q2::AbstractVector,R1::AbstractVector,R2::AbstractVector; qQ=0, qR=0, sysc=ss(0), L=Matrix(), K=Matrix(), integrator=false)
    Q1 = diagm(Q1)
    Q2 = diagm(Q2)
    R1 = diagm(R1)
    R2 = diagm(R2)
    LQG(P,Q1,Q2,R1,R2, qQ, qR, sysc, L, K, integrator)
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
        Acl = [A-B*L B*L; zeros(A) A-K*C]
        Bcl = [B; zeros(B)]
        Ccl = [C zeros(C)]
        Bcl = Bcl/(P.C*inv(P.B*L[:,1:n]-P.A)*P.B) # B*lᵣ # Always normalized with nominal plant static gain
        return syscl = ss(Acl,Bcl,Ccl,0)
    elseif s ∈ [:Sin, :S] # Sensitivity function
        return feedback(ss(eye(m)),PC)
    elseif s ∈ [:Tin, :T] # Complementary sensitivity function
        return feedback(PC)
    elseif s == :Sout # Sensitivity function, output
        return feedback(ss(eye(m)),sysc*P)
    elseif s == :Tout # Complementary sensitivity function, output
        return feedback(sysc*P)
    elseif s == :PS # Load disturbance to output
        return P*G[:S]
    elseif s == :CS # Noise to control signal
        return sysc*G[:S]
    elseif s ∈ [:lt, :looptransfer, :loopgain]
        return PC
    elseif s ∈ [:rd, :returndifference]
        return  ss(eye(p)) + PC
    elseif s ∈ [:sr, :stabilityrobustness]
        return  ss(eye(p)) + inv(PC)
    end
    error("The symbol $s does not have a function associated with it.")
end




function gangoffour(G::LQG)
    return G[:S], G[:PS], G[:CS], G[:T]
end

function gangoffourplot(G::LQG, args...)
    S,D,N,T = gangoffour(G)
    f1 = sigmaplot(S, args...); Plots.plot!(title="\$S = 1/(1+PC)\$")
    f2 = sigmaplot(D, args...); Plots.plot!(title="\$D = P/(1+PC)\$")
    f3 = sigmaplot(N, args...); Plots.plot!(title="\$N = C/(1+PC)\$")
    f4 = sigmaplot(T, args...); Plots.plot!(title="\$T = PC/(1+PC\$)")
    Plots.subplot(f1,f2,f3,f4)
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
