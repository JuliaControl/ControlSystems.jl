export rstd, rstc, dab, c2d_roots2poly, c2d_poly2poly, zpconv#, lsima, indirect_str


"""
    `sysd, x0map = c2d(sys::StateSpace, Ts, method=:zoh)`
    `sysd = c2d(sys::TransferFunction, Ts, method=:zoh)`

Convert the continuous system `sys` into a discrete system with sample time
`Ts`, using the provided method. Currently only `:zoh` and `:foh` are provided.

Returns the discrete system `sysd`, and for StateSpace systems a matrix `x0map` that transforms the
initial conditions to the discrete domain by
`x0_discrete = x0map*[x0; u0]`"""
function c2d(sys::StateSpace, Ts::Real, method::Symbol=:zoh)
    if isdiscrete(sys)
        error("sys must be a continuous time system")
    end
    A, B, C, D = ssdata(sys)
    ny, nu = size(sys)
    nx = nstates(sys)
    if method === :zoh
        M = exp([A*Ts  B*Ts;
            zeros(nu, nx + nu)])
        Ad = M[1:nx, 1:nx]
        Bd = M[1:nx, nx+1:nx+nu]
        Cd = C
        Dd = D
        x0map = [Matrix{Float64}(I, nx, nx) zeros(nx, nu)] # Cant use I if nx==0
    elseif method === :foh
        M = exp([A*Ts B*Ts zeros(nx, nu);
            zeros(nu, nx + nu) Matrix{Float64}(I, nu, nu);
            zeros(nu, nx + 2*nu)])
        M1 = M[1:nx, nx+1:nx+nu]
        M2 = M[1:nx, nx+nu+1:nx+2*nu]
        Ad = M[1:nx, 1:nx]
        Bd = Ad*M2 + M1 - M2
        Cd = C
        Dd = D + C*M2
        x0map = [Matrix{Float64}(I, nx, nx)  (-M2)]
    elseif method === :tustin || method === :matched
        error("NotImplemented: Only `:zoh` and `:foh` implemented so far")
    else
        error("Unsupported method: ", method)
    end
    return StateSpace(Ad, Bd, Cd, Dd, Ts), x0map
end


function d2c(sys::AbstractStateSpace{<:Discrete}, method::Symbol=:zoh)
    A, B, C, D = ssdata(sys)
    ny, nu = size(sys)
    nx = nstates(sys)
    if method == :zoh
        M = log([A  B;
            zeros(nu, nx) I])./sys.Ts
        Ac = M[1:nx, 1:nx]
        Bc = M[1:nx, nx+1:nx+nu]
        Cc = C
        Dc = D
        if eltype(A) <: Real
            Ac,Bc = real.((Ac, Bc))
        end
    else
        error("Unsupported method: ", method)
    end
    return StateSpace(Ac, Bc, Cc, Dc)
end

d2c(sys::TransferFunction{<:Discrete}, args...) = tf(d2c(ss(sys), args...))


function rst(bplus,bminus,a,bm1,am,ao,ar=[1],as=[1] ;cont=true)

    ae      = conv(a,ar)
    be      = conv(bminus,as)
    aoam    = conv(am,ao)
    r1,s1   = dab(ae,be,aoam)

    r       = conv(conv(r1,ar),bplus)
    s       = conv(s1,as)

    bm      = conv(bminus,bm1)
    t0      = (cont ? am[end]/bm[end] : sum(am)/sum(bm))
    t       = t0*conv(ao,bm1)
    s       = s/r[1]
    t       = t/r[1]
    r       = r/r[1]

    r,s,t
end



"""
See ?rstd for the discerte case
"""
rstc(args...)=rst(args..., ;cont=true)

"""
rstd  Polynomial synthesis in discrete time.

`R,S,T=rstd(BPLUS,BMINUS,A,BM1,AM,AO,AR,AS)`

`R,S,T=rstd(BPLUS,BMINUS,A,BM1,AM,AO,AR)`

`R,S,T=rstd(BPLUS,BMINUS,A,BM1,AM,AO)`

Polynomial synthesis according to CCS ch 10 to
design a controller R(q) u(k) = T(q) r(k) - S(q) y(k)

Inputs:  BPLUS  : Part of open loop numerator
BMINUS : Part of open loop numerator
A      : Open loop denominator
BM1    : Additional zeros
AM     : Closed loop denominator
AO     : Observer polynomial
AR     : Pre-specified factor of R,
e.g integral part [1, -1]^k
AS     : Pre-specified factor of S,
e.g notch filter [1, 0, w^2]

Outputs: R,S,T  : Polynomials in controller

See function DAB how the solution to the Diophantine-
Aryabhatta-Bezout identity is chosen.

See Computer-Controlled Systems: Theory and Design, Third Edition
Karl Johan Åström, Björn Wittenmark
"""
rstd(args...)=rst(args..., ;cont=false)


"""
DAB   Solves the Diophantine-Aryabhatta-Bezout identity

`X,Y = DAB(A,B,C)`

AX + BY = C, where A, B, C, X and Y are polynomials
and deg Y = deg A - 1.

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
    length(s) > length(r) && @warn("Controller not casual, deg(S) > deg(R), consider increasing degree of observer polynomial")
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
`c2d_roots2poly(ro,h)`

returns the polynomial coefficients in discrete time given a vector of roots in continuous time
"""
function c2d_roots2poly(ro,h)
    return real((Polynomials.poly(exp(ro.*h))).coeffs[end:-1:1])
end

"""
`c2d_poly2poly(ro,h)`

returns the polynomial coefficients in discrete time given polynomial coefficients in continuous time
"""
function c2d_poly2poly(p,h)
    ro = Polynomials.roots(Polynomials.Polynomial(p[end:-1:1]))
    return real(Polynomials.poly(exp(ro.*h)).coeffs[end:-1:1])
end


function c2d(G::TransferFunction, h, args...)
    @assert iscontinuous(G)
    ny, nu = size(G)
    @assert (ny + nu == 2) "c2d(G::TransferFunction, h) not implemented for MIMO systems"
    sys = ss(G)
    sysd = c2d(sys, h, args...)[1]
    return convert(TransferFunction, sysd)
end

"""
`zpc(a,r,b,s)` form conv(a,r) + conv(b,s) where the lengths of the polynomials are equalized by zero-padding such that the addition can be carried out
"""
function zpconv(a,r,b,s)
    d = length(a)+length(r)-length(b)-length(s)
    if d > 0
        b = [zeros(d);b]
    end
    conv(a,r) + conv(b,s)
end
