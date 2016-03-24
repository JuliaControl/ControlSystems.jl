export rstd, rstc, dab, c2d_roots2poly, c2d_poly2poly, zpconv#, lsima, indirect_str


@doc """`[sysd, x0map] = c2d(sys, Ts, method=:zoh)`

Convert the continuous system `sys` into a discrete system with sample time
`Ts`, using the provided method. Currently only `:zoh` and `:foh` are provided.

Returns the discrete system `sysd`, and a matrix `x0map` that transforms the
initial conditions to the discrete domain by
`x0_discrete = x0map*[x0; u0]`""" ->
function c2d(sys::StateSpace, Ts::Real, method::Symbol=:zoh)
    if !iscontinuous(sys)
        error("sys must be a continuous time system")
    end
    A, B, C, D = sys.A, sys.B, sys.C, sys.D
    ny, nu = size(sys)
    nx = sys.nx
    if method == :zoh
        M = expm([A*Ts  B*Ts;
            zeros(nu, nx + nu)])
        Ad = M[1:nx, 1:nx]
        Bd = M[1:nx, nx+1:nx+nu]
        Cd = C
        Dd = D
        x0map = [eye(nx) zeros(nx, nu)]
    elseif method == :foh
        M = expm([A*Ts B*Ts zeros(nx, nu);
            zeros(nu, nx + nu) eye(nu);
            zeros(nu, nx + 2*nu)])
        M1 = M[1:nx, nx+1:nx+nu]
        M2 = M[1:nx, nx+nu+1:nx+2*nu]
        Ad = M[1:nx, 1:nx]
        Bd = Ad*M2 + M1 - M2
        Cd = C
        Dd = D + C*M2
        x0map = [eye(nx) -M2]
    elseif method == :tustin || method == :matched
        error("NotImplemented: Only `:zoh` and `:foh` implemented so far")
    else
        error("Unsupported method: ", method)
    end
    return StateSpace(Ad, Bd, Cd, Dd, Ts, sys.statenames, sys.inputnames,
        sys.outputnames), x0map
end


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
        warn("Singular problem due to common factors in A and B")
    end
    co = cond(m)
    co > 1e6 && println("dab: condition number $(co)")
    rs = c'/(m')
    r = rs[1:nr]
    s = rs[nr+1:nc]
    length(s) > length(r) && warn("Controller not casual, deg(S) > deg(R), consider increasing degree of observer polynomial")
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

import Polynomials

"""
`c2d_roots2poly(ro,h)`

returns the polynomial coefficients in discrete time given a vector of roots in continuous time
"""
function c2d_roots2poly(ro,h)
    return real((Polynomials.poly(exp(ro.*h))).a[end:-1:1])
end

"""
`c2d_poly2poly(ro,h)`

returns the polynomial coefficients in discrete time given polynomial coefficients in continuous time
"""
function c2d_poly2poly(p,h)
    ro = Polynomials.roots(Polynomials.Poly(p[end:-1:1]))
    return real(Polynomials.poly(exp(ro.*h)).a[end:-1:1])
end


function c2d(G::TransferFunction, h;kwargs...)
    @assert iscontinuous(G)
    ny, nu = size(G)
    @assert (ny + nu == 2) "c2d(G::TransferFunction, h) not implemented for MIMO systems"
    sys = ss(G)
    sysd = c2d(sys,h,kwargs...)[1]
    return ss2tf(sysd)
end


@doc """`[y, t, x] = lsima(sys, t[, x0, method])`

Calculate the time response of adaptive controller. If `x0` is ommitted,
a zero vector is used.

Continuous time systems are discretized before simulation. By default, the
method is chosen based on the smoothness of the input signal. Optionally, the
`method` parameter can be specified as either `:zoh` or `:foh`.""" ->
function lsima{T}(sys::StateSpace, t::AbstractVector, r::AbstractVector{T}, control_signal::Function,state,
    x0::VecOrMat=zeros(sys.nx, 1), method::Symbol=:zoh)
    ny, nu = size(sys)

    nx = sys.nx

    if length(x0) != nx
        error("size(x0) must match the number of states of sys")
    end

    dt = Float64(t[2] - t[1])
    if !iscontinuous(sys) || method == :zoh
        if iscontinuous(sys)
            dsys = c2d(sys, dt, :zoh)[1]
        else
            if sys.Ts != dt
                error("Time vector must match sample time for discrete system")
            end
            dsys = sys
        end
    else
        dsys, x0map = c2d(sys, dt, :foh)
    end
    n = size(t, 1)
    x = Array(T, size(sys.A, 1), n)
    u = Array(T, n)
    y = Array(T, n)
    for i=1:n
        x[:,i] = x0
        y[i] = (sys.C*x0 + sys.D*u[i])[1]

        u[i],state = control_signal(state, y[1:i], u[1:i-1], r[1:i])
        x0 = sys.A * x0 + sys.B * u[i]

    end


    return y, t, x, u
end
lsima(sys::TransferFunction, u, t,r, args...) = lsima(ss(sys), u, t,r, args...)

function indirect_str(state, y, u,uc, nb,na, lambda,bm1,am,ao,ar=[1],as=[1])
    theta, P = state
    u = [zeros(length(ao)+length(am)+1);u]
    y = [zeros(length(ao)+length(am)+1);y]
    phi = [y[end-1:-1:end-na]; u[end:-1:end-nb]]
    # compute new estimate and update covariance matrix

    K = P*phi/(lambda + phi'P*phi)
    new_theta = theta + K*(y[end] - phi'theta)
    new_P = (I - K*phi')*P/lambda
    new_P = (new_P + new_P')/2

    state = (new_theta, new_P)

    a = [1;theta[1:na]]
    b = theta[na+1:end]
    r,s,t = rstd([1],b,a,bm1,am,ao,ar,as)
    uo = r⋅u[end:-1:end-length(r)+1] + s⋅y[end-1:-1:end-length(s)] + t⋅uc[end:-1:end-length(t)+1]

    return uo,state

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
