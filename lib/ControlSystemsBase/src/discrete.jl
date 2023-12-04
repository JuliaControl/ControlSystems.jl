export rstd, rstc, dab, c2d_roots2poly, c2d_poly2poly, zpconv#, lsima, indirect_str
using LinearAlgebra: exp!


"""
    sysd = c2d(sys::AbstractStateSpace{<:Continuous}, Ts, method=:zoh; w_prewarp=0)
    Gd = c2d(G::TransferFunction{<:Continuous}, Ts, method=:zoh)

Convert the continuous-time system `sys` into a discrete-time system with sample time
`Ts`, using the specified `method` (:`zoh`, `:foh`, `:fwdeuler` or `:tustin`).

`method = :tustin` performs a bilinear transform with prewarp frequency `w_prewarp`.

- `w_prewarp`: Frequency (rad/s) for pre-warping when using the Tustin method, has no effect for other methods.

See also `c2d_x0map`

# Extended help

ZoH sampling is exact for linear systems with piece-wise constant inputs (step invariant), i.e., the solution obtained using [`lsim`](@ref) is not approximative (modulu machine precision). ZoH sampling is commonly used to discretize continuous-time plant models that are to be controlled using a discrete-time controller.

FoH sampling is exact for linear systems with piece-wise linear inputs (ramp invariant), this is a good choice for simulation of systems with smooth continuous inputs.

To approximate the behavior of a continuous-time system well in the frequency domain, the `:tustin` (trapezoidal / bilinear) method may be most appropriate. In this case, the pre-warping argument can be used to ensure that the frequency response of the discrete-time system matches the continuous-time system at a given frequency. The tustin transformation alters the meaning of the state components, while ZoH and FoH preserve the meaning of the state components. The Tustin method is commonly used to discretize a continuous-tiem controller.

The forward-Euler method generally requires the sample time to be very small
relative to the time constants of the system, and its use is generally discouraged.

Classical rules-of-thumb for selecting the sample time for control design dictate
that `Ts` should be chosen as ``0.2 ≤ ωgc⋅Ts ≤ 0.6`` where ``ωgc`` is the gain-crossover frequency (rad/s).
"""
c2d(sys::AbstractStateSpace{<:Continuous}, Ts::Real, method::Symbol=:zoh; kwargs...) = c2d_x0map(sys, Ts, method; kwargs...)[1]


"""
    sysd, x0map = c2d_x0map(sys::AbstractStateSpace{<:Continuous}, Ts, method=:zoh; w_prewarp=0)

Returns the discretization `sysd` of the system `sys` and a matrix `x0map` that
transforms the initial conditions to the discrete domain by `x0_discrete = x0map*[x0; u0]`

See `c2d` for further details.
"""
function c2d_x0map(sys::AbstractStateSpace{<:Continuous}, Ts::Real, method::Symbol=:zoh; w_prewarp=0)
    A, B, C, D = ssdata(sys)
    T = promote_type(eltype.((A,B,C,D))...)
    ny, nu = size(sys)
    nx = nstates(sys)
    if method === :zoh
        M = exp!([A*Ts  B*Ts;
            zeros(nu, nx + nu)])
        Ad = M[1:nx, 1:nx]
        Bd = M[1:nx, nx+1:nx+nu]
        Cd = C
        Dd = D
        x0map = [Matrix{T}(I, nx, nx) zeros(nx, nu)] # Can't use I if nx==0
    elseif method === :foh
        M = exp!([A*Ts B*Ts zeros(nx, nu);
            zeros(nu, nx + nu) Matrix{T}(I, nu, nu);
            zeros(nu, nx + 2*nu)])
        M1 = M[1:nx, nx+1:nx+nu]
        M2 = M[1:nx, nx+nu+1:nx+2*nu]
        Ad = M[1:nx, 1:nx]
        Bd = Ad*M2 + M1 - M2
        Cd = C
        Dd = D + C*M2
        x0map = [Matrix{T}(I, nx, nx)  (-M2)]
    elseif method === :fwdeuler
        Ad, Bd, Cd, Dd = (I+Ts*A), Ts*B, C, D
        x0map = I(nx)
    elseif method === :tustin
        a = w_prewarp == 0 ? Ts/2 : tan(w_prewarp*Ts/2)/w_prewarp
        a > 0 || throw(DomainError("A positive w_prewarp must be provided for method Tustin"))
        AI = (I-a*A)
        Ad = AI\(I+a*A)
        Bd = 2a*(AI\B)
        Cd = C/AI 
        Dd = a*Cd*B + D
        x0map = Matrix{T}(I, nx, nx)
    elseif method === :matched
        error("NotImplemented: Only `:zoh`, `:foh`, :tustin and `:fwdeuler` implemented so far")
    else
        error("Unsupported method: ", method)
    end
    timeevol = Discrete(Ts)
    return StateSpace{typeof(timeevol), eltype(Ad)}(Ad, Bd, Cd, Dd, timeevol), x0map
end

"""
    d2c(sys::AbstractStateSpace{<:Discrete}, method::Symbol = :zoh; w_prewarp=0)

Convert discrete-time system to a continuous time system, assuming that the discrete-time system was discretized using `method`. Available methods are `:zoh, :fwdeuler´.

- `w_prewarp`: Frequency for pre-warping when using the Tustin method, has no effect for other methods.

See also [`d2c_exact`](@ref).
"""
function d2c(sys::AbstractStateSpace{<:Discrete}, method::Symbol=:zoh; w_prewarp=0)
    A, B, C, D = ssdata(sys)
    ny, nu = size(sys)
    nx = nstates(sys)
    if method === :zoh
        M = log([A  B;
            zeros(nu, nx) I])./sys.Ts
        Ac = M[1:nx, 1:nx]
        Bc = M[1:nx, nx+1:nx+nu]
        if eltype(A) <: Real
            Ac,Bc = real.((Ac, Bc))
        end
        Cc, Dc = C, D
    elseif method === :fwdeuler
        Ac = (A-I)./sys.Ts
        Bc = B./sys.Ts
        Cc, Dc = C, D
    elseif method === :tustin
        a = w_prewarp == 0 ? sys.Ts/2 : tan(w_prewarp*sys.Ts/2)/w_prewarp
        a > 0 || throw(DomainError("A positive w_prewarp must be provided for method Tustin"))
        AI = a*(A+I)
        Ac = (A-I)/AI
        Bc = AI\B
        Cc = 2a*C/AI
        Dc = D - Cc*B/2
    else
        error("Unsupported method: ", method)
    end
    return StateSpace(Ac, Bc, Cc, Dc)
end

d2c(sys::TransferFunction{<:Discrete}, args...) = tf(d2c(ss(sys), args...))


"""
    d2c_exact(sys::AbstractStateSpace{<:Discrete})

Translate a discrete-time system to a continuous-time system by the substitution ``z = e^{sT_s}``.
The translation is exact in the frequency domain, i.e.,
the frequency response of the resulting continuous-time system is identical to
the frequency response of the discrete-time system.

This method is useful when analyzing hybrid continuous/discrete systems in the frequency domain and high accuracy is required.

The resulting system will be be a static system in feedback with pure negative delays,
i.e., this system cannot be simulated in the time domain.
"""
function d2c_exact(sys::AbstractStateSpace{<:Discrete})
    T = sys.Ts
    A,B,C,D = ssdata(sys)
    z = delay(-T)
    LR = append([z for _ in 1:sys.nx]...) - ss(A + I)
    C*feedback(I(sys.nx), LR)*B + D
end

# c2d and d2c for covariance and cost matrices =================================
"""
    Qd     = c2d(sys::StateSpace{Continuous}, Qc::Matrix, Ts;             opt=:o)
    Qd, Rd = c2d(sys::StateSpace{Continuous}, Qc::Matrix, Rc::Matrix, Ts; opt=:o)
    Qd     = c2d(sys::StateSpace{Discrete},   Qc::Matrix;                 opt=:o)
    Qd, Rd = c2d(sys::StateSpace{Discrete},   Qc::Matrix, Rc::Matrix;     opt=:o)

Sample a continuous-time covariance or LQR cost matrix to fit the provided discrete-time system.

If `opt = :o` (default), the matrix is assumed to be a covariance matrix. The measurement covariance `R` may also be provided.
If `opt = :c`, the matrix is instead assumed to be a cost matrix for an LQR problem.

!!! note
    Measurement covariance (here called `Rc`) is usually estimated in discrete time, and is in this case not dependent on the sample rate. Discretization of the measurement covariance only makes sense when a continuous-time controller has been designed and the closest corresponding discrete-time controller is desired.

The method used comes from theorem 5 in the reference below.

Ref: "Discrete-time Solutions to the Continuous-time
Differential Lyapunov Equation With Applications to Kalman Filtering", 
Patrik Axelsson and Fredrik Gustafsson

On singular covariance matrices: The traditional double integrator with covariance matrix `Q = diagm([0,σ²])` can not be sampled with this method. Instead, the input matrix ("Cholesky factor") of `Q` must be manually kept track of, e.g., the noise of variance `σ²` enters like `N = [0, 1]` which is sampled using ZoH and becomes `Nd = [1/2 Ts^2; Ts]` which results in the covariance matrix `σ² * Nd * Nd'`. 

# Example:
The following example designs a continuous-time LQR controller for a resonant system. This is simulated with OrdinaryDiffEq to allow the ODE integrator to also integrate the continuous-time LQR cost (the cost is added as an additional state variable). We then discretize both the system and the cost matrices and simulate the same thing. The discretization of an LQR contorller in this way is sometimes refered to as `lqrd`.
```julia
using ControlSystemsBase, LinearAlgebra, OrdinaryDiffEq, Test
sysc = DemoSystems.resonant()
x0 = ones(sysc.nx)
Qc = [1 0.01; 0.01 2] # Continuous-time cost matrix for the state
Rc = I(1)             # Continuous-time cost matrix for the input

L = lqr(sysc, Qc, Rc)
dynamics = function (xc, p, t)
    x = xc[1:sysc.nx]
    u = -L*x
    dx = sysc.A*x + sysc.B*u
    dc = dot(x, Qc, x) + dot(u, Rc, u)
    return [dx; dc]
end
prob = ODEProblem(dynamics, [x0; 0], (0.0, 10.0))
sol = solve(prob, Tsit5(), reltol=1e-8, abstol=1e-8)
cc = sol.u[end][end] # Continuous-time cost

# Discrete-time version
Ts = 0.01 
sysd = c2d(sysc, Ts)
Ld = lqr(sysd, Qd, Rd)
sold = lsim(sysd, (x, t) -> -Ld*x, 0:Ts:10, x0 = x0)
function cost(x, u, Q, R)
    dot(x, Q, x) + dot(u, R, u)
end
cd = cost(sold.x, sold.u, Qd, Rd) # Discrete-time cost
@test cc ≈ cd rtol=0.01           # These should be similar
```
"""
function c2d(sys::AbstractStateSpace{<:Continuous}, Qc::AbstractMatrix, Ts::Real; opt=:o)
    n = sys.nx
    Ac  = sys.A
    if opt === :c
        # Ref: Charles Van Loan: Computing integrals involving the matrix exponential, IEEE Transactions on Automatic Control. 23 (3): 395–404, 1978
        F = [-Ac' Qc; zeros(size(Qc)) Ac]
        G = exp(F*Ts)
        Ad = G[n+1:end, n+1:end]'
        AdiQd = G[1:n, n+1:end]
        Qd = Ad*AdiQd
    elseif opt === :o
        # Ref: Discrete-time Solutions to the Continuous-time Differential Lyapunov Equation With Applications to Kalman Filtering
        F = [Ac Qc; zeros(size(Qc)) -Ac']
        M = exp(F*Ts)
        M1 = M[1:n, 1:n]
        M2 = M[1:n, n+1:end]
        Qd = M2*M1'
    else
        error("Unknown option opt=$opt")
    end
    
    (Qd .+ Qd') ./ 2
end


function c2d(sys::AbstractStateSpace{<:Discrete}, Qc::AbstractMatrix, R::Union{AbstractMatrix, Nothing}=nothing; opt=:o)
    Ad  = sys.A
    Ac  = real(log(Ad)./sys.Ts)
    if opt === :c
        Ac = Ac'
        Ad = Ad'
    elseif opt !== :o
        error("Unknown option opt=$opt")
    end
    C   = Symmetric(Qc - Ad*Qc*Ad')
    Qd  = MatrixEquations.lyapc(Ac, C)
    # The method below also works, but no need to use quadgk when MatrixEquations is available.
    # function integrand(t)
    #     Ad = exp(t*Ac)
    #     Ad*Qc*Ad'
    # end
    # Qd = quadgk(integrand, 0, h)[1]
    if R === nothing
        return Qd
    else
        if opt === :c
            Qd, R .* sys.Ts
        else
            Qd, R ./ sys.Ts
        end
    end
end


function c2d(sys::AbstractStateSpace, Qc::AbstractMatrix, R::AbstractMatrix, Ts::Real; opt=:o)
    Qd = c2d(sys, Qc, Ts; opt)
    if opt === :c
        return Qd, R .* Ts
    else
        return Qd, R ./ Ts
    end
end



"""
    Qc = d2c(sys::AbstractStateSpace{<:Discrete}, Qd::AbstractMatrix; opt=:o)

Resample discrete-time covariance matrix belonging to `sys` to the equivalent continuous-time matrix.

The method used comes from theorem 5 in the reference below.

If `opt = :c`, the matrix is instead assumed to be a cost matrix for an LQR problem.

Ref: Discrete-time Solutions to the Continuous-time
Differential Lyapunov Equation With
Applications to Kalman Filtering
Patrik Axelsson and Fredrik Gustafsson
"""
function d2c(sys::AbstractStateSpace{<:Discrete}, Qd::AbstractMatrix, Rd::Union{AbstractMatrix, Nothing}=nothing; opt=:o)
    Ad = sys.A
    Ac = real(log(Ad)./sys.Ts)
    if opt === :c
        Ac = Ac'
        Ad = Ad'
    elseif opt !== :o
        error("Unknown option opt=$opt")
    end
    C = Symmetric(Ac*Qd + Qd*Ac')
    Qc = MatrixEquations.lyapd(Ad, -C)
    isposdef(Qc) || @error("Calculated covariance matrix not positive definite")
    if Rd === nothing
        return Qc
    else
        if opt === :c
            return Qc, Rd ./ sys.Ts
        else
            return Qc, Rd .* sys.Ts
        end
    end
end


function rst(bplus,bminus,a,bm1,am,ao,ar=[1],as=[1] ;cont=true)

    ae      = conv(a,ar)
    be      = conv(bminus,as)
    aoam    = conv(am,ao)
    ret   = dab(ae,be,aoam)
    ret === nothing && error("Controller not casual, deg(S) > deg(R), consider increasing degree of observer polynomial")
    r1,s1 = ret

    r       = conv(conv(r1,ar),bplus)
    s       = conv(s1,as)

    bm      = conv(bminus,bm1)
    t0      = (cont ? am[end]/bm[end] : sum(am)/sum(bm))
    t       = t0*conv(ao,bm1)
    s       = s/r[1]
    t       = t/r[1]
    r       = r/r[1]

    length(s) > length(r) && @warn("Controller not casual, deg(S) > deg(R), consider increasing degree of observer polynomial")

    r,s,t
end



"""
See `?rstd` for the discrete case
"""
rstc(args...)=rst(args..., ;cont=true)

"""
    R,S,T = rstd(BPLUS,BMINUS,A,BM1,AM,AO,AR,AS)
    R,S,T = rstd(BPLUS,BMINUS,A,BM1,AM,AO,AR)
    R,S,T = rstd(BPLUS,BMINUS,A,BM1,AM,AO)

Polynomial synthesis in discrete time.

Polynomial synthesis according to "Computer-Controlled Systems" ch 10 to
design a controller ``R(q) u(k) = T(q) r(k) - S(q) y(k)``

Inputs:
- `BPLUS`  : Part of open loop numerator
- `BMINUS` : Part of open loop numerator
- `A`      : Open loop denominator
- `BM1`    : Additional zeros
- `AM`     : Closed loop denominator
- `AO`     : Observer polynomial
- `AR`     : Pre-specified factor of R,
e.g integral part [1, -1]^k
- `AS`     : Pre-specified factor of S,
e.g notch filter [1, 0, w^2]

Outputs: `R,S,T`  : Polynomials in controller

See function `dab` how the solution to the Diophantine-
Aryabhatta-Bezout identity is chosen.

See Computer-Controlled Systems: Theory and Design, Third Edition
Karl Johan Åström, Björn Wittenmark
"""
rstd(args...)=rst(args..., ;cont=false)


"""
    X,Y = dab(A,B,C)

Solves the Diophantine-Aryabhatta-Bezout identity

``AX + BY = C``, where ``A, B, C, X`` and ``Y`` are polynomials
and ``deg Y = deg A - 1``.

See Computer-Controlled Systems: Theory and Design, Third Edition
Karl Johan Åström, Björn Wittenmark
"""
function dab(a,b,c)

    na = length(a)
    nb = length(b)
    nc = length(c)
    ns = na - 1
    if ns < 1
        r = c/a
        s = 0
        return
    end
    nr = nc - ns
    c = nb-nr > 1 ? [zeros(nb-nr-1); c] : c
    nc = length(c)
    nr = nc - ns
    if nr < 1
        r = 0
        s = c/b
        return
    end
    b = [zeros(nr-nb+1); b]
    za = zeros(nr-1)
    zb = zeros(ns-1)
    ma = toeplitz([a; za],[a[1]; za])
    mb = toeplitz([b; zb],[b[1]; zb])
    m = [ma mb]
    if rank(m) < minimum(size(m))
        @warn("Singular problem due to common factors in A and B")
    end
    co = cond(m)
    co > 1e6 && println("dab: condition number $(co)")
    rs = (c'/(m'))'
    r = rs[1:nr]
    s = rs[nr+1:nc]
    r,s
end

function toeplitz(c,r)
    nc = length(c)
    nr = length(r)
    A  = zeros(nc, nr)
    A[:,1] = c
    A[1,:] = r
    for i in 2:nr
        A[2:end,i] = A[1:end-1,i-1]
    end
    A
end


"""
    c2d_roots2poly(ro, Ts)

returns the polynomial coefficients in discrete time given a vector of roots in continuous time
"""
function c2d_roots2poly(ro, Ts)
    return reverse(Polynomials.fromroots(exp.(ro .* Ts)).coeffs)
end

"""
    c2d_poly2poly(ro, Ts)

returns the polynomial coefficients in discrete time given polynomial coefficients in continuous time
"""
function c2d_poly2poly(p, Ts)
    ro = Polynomials.roots(Polynomials.Polynomial(reverse(p)))
    return c2d_roots2poly(ro, Ts)
end

function c2d(G::TransferFunction{<:Continuous, <:SisoRational}, Ts, args...; kwargs...)
    issiso(G) || error("c2d(G::TransferFunction, h) not implemented for MIMO systems")
    sysd = c2d(ss(G), Ts, args...; kwargs...)
    return convert(TransferFunction{typeof(sysd.timeevol), SisoRational}, sysd)
end

function c2d(G::TransferFunction{<:Continuous, <:SisoZpk}, Ts, args...; kwargs...)
    issiso(G) || error("c2d(G::TransferFunction, h) not implemented for MIMO systems")
    sysd = c2d(ss(G), Ts, args...; kwargs...)
    return convert(TransferFunction{typeof(sysd.timeevol), SisoZpk}, sysd)
end

"""
    zpc(a,r,b,s)

form `conv(a,r) + conv(b,s)` where the lengths of the polynomials are equalized by zero-padding such that the addition can be carried out
"""
function zpconv(a,r,b,s)
    d = length(a)+length(r)-length(b)-length(s)
    if d > 0
        b = [zeros(d);b]
    end
    conv(a,r) + conv(b,s)
end
