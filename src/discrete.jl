@doc """`[sysd, x0map] = c2d(sys, Ts, method=:zoh)`

Convert the continuous system `sys` into a discrete system with sample time
`Ts`, using the provided method. Currently only `:zoh` and `:foh` are provided.

Returns the discrete system `sysd`, and a matrix `x0map` that transforms the
initial conditions to the discrete domain by
`x0_discrete = x0map*[x0; u0]`"""->
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
