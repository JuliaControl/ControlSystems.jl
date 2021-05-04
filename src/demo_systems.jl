"""
    ssrand(T::Type, ny::Int, nu::Int, nstates::Int; proper=false, stable=true, Ts=nothing)

Returns a random `StateSpace` model with `ny` outputs, `nu` inputs, and `nstates` states,
whose matrix elements are normally distributed.

It is possible to specify if the system should be `proper` or `stable`.

Specify a sample time `Ts` to obtain a discrete-time system.
"""
function ssrand(T::Type, ny::Int, nu::Int, nstates::Int; proper=false, stable=true, Ts=nothing)
    A = randn(T, nstates, nstates)
    if stable
        Λ = eigvals(A)
        A = isnothing(Ts) ? A - 1.1*max(0, maximum(real(Λ)))*I : A*0.9/maximum(abs.(Λ))
    end

    B = randn(T, nstates, nu)
    C = randn(T, ny, nstates)
    D = (proper ? zeros(T, ny, nu) : randn(T, ny, nu))
    return isnothing(Ts) ? StateSpace(A, B, C, D) : StateSpace(A, B, C, D, Ts)
end
ssrand(ny::Int, nu::Int, nstates::Int; kwargs...) = ssrand(Float64, ny::Int, nu::Int, nstates::Int; kwargs...)
ssrand(n; kwargs...) = ssrand(n, n, 2*n; kwargs...)
ssrand(T::Type, n; kwargs...) = ssrand(T, n, n, 2*n; kwargs...)


"""
    DemoSystems

A module with standard test systems from the control literature.

The returned systems are of type `StateSpace` or `DelayLTISystem`.

Default parameter can be changed using keyword arguments.

SISO systems
============
* `DemoSystems.lag(;T=1)`     First-order system  `1/(sT+1)`
* `DemoSystems.fotd(;T=1, τ=1)`   First-order system with time delay `exp(-sτ)/(sT+1)`
* `DemoSystems.sotd(;T=1, T2=10, τ=1)`   Second-order non-resonant system with time delay `exp(-sτ)/(sT+1)/(sT2 + 1)`
* `DemoSystems.resonant(;ω0=1, ζ=0.25)`   Second-order resonant systems `ω0^2/(s^2 + 2ζ*ω0*s + ω0^2)`

MIMO systems
============
* `DemoSystems.woodberry()`     Wood--Berry distillation column
* `DemoSystems.doylesat(;a=10)`    The spinning body example by Doyle

Fore more information, see the help for the specific systems.
"""
module DemoSystems

using ControlSystems
using LinearAlgebra


"""
    lag(;T=1, τ=1)

Returns a first-order system `StateSpace` with transfer function `1/(sT+1)`.
"""
lag(;T=1) = ss(-1/T, 1, 1/T, 0)


"""
    fotd(;T=1, τ=1)

Returns a first-order system + a time delay `exp(-sτ)/(sT+1)` (`DelayLTISystem`).
"""
fotd(;T=1, τ=1) = ss(-1/T, 1, 1/T, 0) * delay(τ)


"""
    sotd(;T=1, T2=10, τ=1)

Returns a second-order system + a time delay `exp(-sτ)/(sT+1)/(sT2 + 1)` (`DelayLTISystem`).
"""
sotd(;T=1, T2=10, τ=1) = ss(-1/T, 1, 1/T, 0)*ss(-1/T2, 1, 1/T2, 0)*delay(τ)


"""
    resonant(;ω0=1, ζ=0.25)

Returns a second-order system + a time delay `exp(-sτ)/(sT+1)/(sT2 + 1)` (`DelayLTISystem`).
"""
resonant(;ω0=1, ζ=0.25) = ss([-ζ -ω0; ω0 -ζ], [ω0; 0], [0 ω0], 0) # QUESTION: Is this the form that we like? Perhhaps not.


"""
    `doylesat(;a=10) = ss([0 a; -a 0], I(2), [1 a; -a 1], 0)`

Returns a `StateSpace` model for the spinning body (satellite) example by Doyle,
which is a classic example that illustrates that robustness analysis of MIMO systems
is more involved than for SISO systems.

*References:*

**Zhou, K., J. C. Doyle, and K. Glover**, Robust and optimal control,
Prentice hall (NJ), 1996.

**Skogestad, S, and I. Postlethwaite**, Multivariable feedback control:
analysis and design, Wiley (NY), 2007.
"""
doylesat(;a=10) = ss([0 a; -a 0], I(2), [1 a; -a 1], 0)


"""
    `woodberry()`

Returns a `DelayLTISystem` model for the Wood--Berry distillation column,
which is a classic example from the literature on process control of MIMO process.

*References:*

**Wood, R. K., and M. W. Berry**, Terminal composition control of a binary
distillation column, Chemical Engineering Science, 28.9, 1973, pp. 1707-1717.

"""
function woodberry()
    s = tf("s")
    return [12.8/(1 + 16.7*s)*delay(1.0)   -18.9/(1 + 21*s) * delay(3.0)
            6.6/(1 + 10.9*s) * delay(7.0)  -19.4/(1 + 14.4*s) * delay(3.0)]
end


end
