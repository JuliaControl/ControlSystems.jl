struct Saturation{T} <: Function
    l::T
    u::T
end
Saturation(u) = Saturation(-u, u)

(s::Saturation)(x) = clamp(x, s.l, s.u)

"""
    saturation(val)
    saturation(lower, upper)

Create a saturating nonlinearity. Connect it to the output of a controller `C` using
```
Csat = saturation(val) * C
```

```
           y▲   ────── upper
            │  /
            │ /
            │/
  ──────────┼────────► u
           /│   
          / │
         /  │
lower──── 
```

$nonlinear_warning

Note: when composing linear systems with nonlinearities, it's often important to handle operating points correctly.
See [`ControlSystemsBase.offset`](@ref) for handling operating points.
"""
saturation(args...) = nonlinearity(Saturation(args...))
saturation(v::AbstractVector, args...) = nonlinearity(Saturation.(v, args...))
Base.show(io::IO, f::Saturation) = f.u == -f.l ? print(io, "saturation($(f.u))") : print(io, "saturation($(f.l), $(f.u))")

## Offset ======================================================================

struct Offset{T} <: Function
    o::T
end

(s::Offset)(x) = x .+ s.o

"""
    offset(val)

Create a constant-offset nonlinearity `x -> x + val`.

$nonlinear_warning

# Example:
To create a linear system that operates around operating point `y₀, u₀`, use
```julia
offset_sys = offset(y₀) * sys * offset(-u₀)
```
note the sign on the offset `u₀`. This ensures that `sys` operates in the coordinates
`Δu = u-u₀, Δy = y-y₀` and the inputs and outputs to the offset system are in their non-offset coordinate system. If the system is linearized around `x₀`, `y₀` is given by `C*x₀`. Additional information and an example is available here https://juliacontrol.github.io/ControlSystemsBase.jl/latest/lib/nonlinear/#Non-zero-operating-point
"""
offset(val::Number) = nonlinearity(Offset(val))
offset(v::AbstractVector) = nonlinearity(Offset.(v))

Base.show(io::IO, f::Offset) = print(io, "offset($(f.o))")



## Ratelimit ===================================================================

"""
    ratelimit(val; Tf)
    ratelimit(lower, upper; Tf)

Create a nonlinearity that limits the rate of change of a signal, roughly equivalent to `` 1/s ∘ sat ∘ s``. `Tf` controls the filter time constant on the derivative used to calculate the rate. 
$nonlinear_warning
"""
function ratelimit(args...; Tf)
    T = promote_type(typeof(Tf), typeof.(args)...)
    # Filtered derivative
    A = [0 1; -2/Tf^2 -2/Tf]
    B = [0; 1]
    C = [4/Tf^2 0]
    F = ss(A,B,C,0)

    integrator = tf([T(Tf), one(T)], [one(T), zero(T)])
    integrator*saturation(args...)*F
end



## DeadZone ====================================================================

struct DeadZone{T} <: Function
    l::T
    u::T
end
DeadZone(u) = DeadZone(-u, u)

function (s::DeadZone)(x)
    if x > s.u
        x - s.u
    elseif x < s.l
        x - s.l
    else
        zero(x)
    end
end

"""
    deadzone(val)
    deadzone(lower, upper)

Create a dead-zone nonlinearity.

```
       y▲
        │     /
        │    /
  lower │   /
─────|──┼──|───────► u
    /   │   upper
   /    │
  /     │
```

$nonlinear_warning

Note: when composing linear systems with nonlinearities, it's often important to handle operating points correctly.
See [`ControlSystemsBase.offset`](@ref) for handling operating points.
"""
deadzone(args...) = nonlinearity(DeadZone(args...))
deadzone(v::AbstractVector, args...) = nonlinearity(DeadZone.(v, args...))
Base.show(io::IO, f::DeadZone) = f.u == -f.l ? print(io, "deadzone($(f.u))") : print(io, "deadzone($(f.l), $(f.u))")


## Hysteresis ==================================================================

"""
    Hysteresis(amplitude, width, hardness)

A named struct representing the internal nonlinearity used inside [`hysteresis`](@ref).
Stores the `amplitude`, `width`, and `hardness` parameters so that analytical
[`describing_function`](@ref) dispatch is possible.

Note: as a callable, `Hysteresis` evaluates the *static* internal function
`width * tanh(hardness * y)`, while the analytical [`describing_function`](@ref)
overload describes the complete relay-with-hysteresis element realized by
[`hysteresis`](@ref) (which embeds this static function in a feedback loop).
"""
struct Hysteresis{T} <: Function
    amplitude::T
    width::T
    hardness::T
end
Hysteresis(amplitude, width, hardness) = Hysteresis(promote(amplitude, width, hardness)...)

function (h::Hysteresis)(y)
    if isfinite(h.hardness)
        h.width * tanh(h.hardness * y)
    else
        h.width * sign(y)
    end
end

Base.show(io::IO, h::Hysteresis) = print(io, "Hysteresis(amplitude=$(h.amplitude), width=$(h.width), hardness=$(h.hardness))")

"""
    hysteresis(; amplitude, width, Tf, hardness)

Create a hysteresis nonlinearity. The signal switches between `±amplitude` when the input crosses `±width`. `Tf` controls the time constant of the internal state that tracks the hysteresis, and `hardness` controls how sharp the transition is between the two states.

```

      y▲
       │
amp┌───┼───┬─►
   │   │   ▲
   │   │   │
 ──┼───┼───┼───► u
   │   │   │w
   ▼   │   │
 ◄─┴───┼───┘
       │
```

$nonlinear_warning
"""
function hysteresis(; amplitude=1.0, width=1.0, Tf=0.001, hardness=20.0)
    T = promote_type(typeof(amplitude), typeof(width), typeof(Tf), typeof(hardness))
    G = tf(1, [Tf, 1])
    nl_func = Hysteresis(T(amplitude), T(width), T(hardness))
    nl = nonlinearity(nl_func)
    amplitude/width*(feedback(G, -nl) - 1)
end


## Describing Function =========================================================

"""
    describing_function(f, A; N=1000)

Compute the describing function ``N(A)`` for a static nonlinearity `f` at input
amplitude `A`. The describing function is the complex ratio of the fundamental
Fourier component of the output to the input amplitude.

The generic method uses numerical integration (trapezoidal rule with `N` points).
The trapezoidal rule converges very rapidly for smooth periodic integrands, but only
slowly when `f` is discontinuous (e.g., a relay); increase `N` if high accuracy is
required for such nonlinearities.

Analytical overloads are provided for [`Saturation`](@ref ControlSystemsBase.saturation),
[`DeadZone`](@ref ControlSystemsBase.deadzone), and
[`Hysteresis`](@ref ControlSystemsBase.Hysteresis).

# Arguments
- `f`: A scalar nonlinear function ``f: \\mathbb{R} \\to \\mathbb{R}``, or a
  [`HammersteinWienerSystem`](@ref) with a single nonlinearity.
- `A`: Input sinusoidal amplitude (positive real).
- `N`: Number of quadrature points for the numerical method.

# Examples
```julia
julia> using ControlSystemsBase

julia> describing_function(Saturation(1.0), 2.0) # Analytical
0.6089977810442294 + 0.0im

julia> describing_function(x -> clamp(x, -1, 1), 2.0) # Numerical
0.6089971740034709 + 2.263467191454538e-17im
```
"""
describing_function(f, A; kwargs...) = describing_function_numerical(f, A; kwargs...)

# Numerical quadrature core. Call this directly instead of `describing_function`
# where analytical dispatch on the nonlinearity type must be bypassed, e.g., to
# compute the describing function of the static callable of a `Hysteresis`.
function describing_function_numerical(f, A; N=1000)
    A > 0 || throw(ArgumentError("Amplitude A must be positive, got $A"))
    dθ = 2π / N
    b1 = 0.0
    a1 = 0.0
    for i in 0:N-1
        θ = i * dθ
        sinθ = sin(θ)
        cosθ = cos(θ)
        val = f(A * sinθ)
        b1 += val * sinθ
        a1 += val * cosθ
    end
    Complex(b1, a1) * dθ / (π * A)
end

"""
    describing_function(s::Saturation, A)

Analytical describing function for symmetric saturation with limit `d`.
For `A ≤ d`, returns `1`. For `A > d`:
```math
N(A) = \\frac{2}{\\pi}\\left[\\arcsin\\frac{d}{A} + \\frac{d}{A}\\sqrt{1 - \\frac{d^2}{A^2}}\\right]
```
Falls back to numerical computation for asymmetric saturation (`l ≠ -u`).
"""
function describing_function(s::Saturation, A; kwargs...)
    A > 0 || throw(ArgumentError("Amplitude A must be positive, got $A"))
    if s.l != -s.u
        return describing_function_numerical(s, A; kwargs...)
    end
    d = s.u
    A ≤ d && return complex(float(one(d)))
    r = d / A
    complex((2/π) * (asin(r) + r * sqrt(1 - r^2)))
end

"""
    describing_function(dz::DeadZone, A)

Analytical describing function for dead-zone. For symmetric dead-zone with threshold
`d` and `A ≤ d`, returns `0`. For `A > d`:
```math
N(A) = 1 - \\frac{2}{\\pi}\\left[\\arcsin\\frac{d}{A} + \\frac{d}{A}\\sqrt{1 - \\frac{d^2}{A^2}}\\right]
```
Computed through the identity `dz(x) = x - clamp(x, l, u)`, which by linearity of the
describing function gives ``N_{dz}(A) = 1 - N_{sat}(A)``. This holds for asymmetric
dead-zone (`l ≠ -u`) as well, in which case the saturation part is computed numerically.
"""
function describing_function(dz::DeadZone, A; kwargs...)
    1 - describing_function(Saturation(dz.l, dz.u), A; kwargs...)
end

"""
    describing_function(h::Hysteresis, A)

Analytical describing function for an ideal relay with hysteresis, i.e., the complete
hysteresis element realized by [`hysteresis`](@ref) (`hardness` is ignored).
For `A ≤ width`, returns `0`. For `A > width`:
```math
N(A) = \\frac{4M}{\\pi A}\\left[\\sqrt{1 - \\frac{b^2}{A^2}} - j\\frac{b}{A}\\right]
```
where `M = amplitude` and `b = width`.

Note: this is *not* the describing function of the static callable
`h(y) = width * tanh(hardness * y)`; to compute that, bypass analytical dispatch
with a closure: `describing_function(y -> h(y), A)`.
"""
function describing_function(h::Hysteresis, A; kwargs...)
    A > 0 || throw(ArgumentError("Amplitude A must be positive, got $A"))
    M = h.amplitude
    b = h.width
    A ≤ b && return complex(zero(float(M)))
    r = b / A
    (4M / (π * A)) * complex(sqrt(1 - r^2), -r)
end

"""
    describing_function(sys::HammersteinWienerSystem, A; kwargs...)

Compute the describing function for the nonlinearity contained in a
[`HammersteinWienerSystem`](@ref). The system must contain exactly one nonlinearity.

The result is the describing function of the static function actually applied by the
system. In particular, for `nonlinearity(Hysteresis(...))` this is the numerically
computed describing function of the static callable `width*tanh(hardness*y)`, not the
analytical relay-with-hysteresis formula of `describing_function(::Hysteresis, A)`.
"""
function describing_function(sys::HammersteinWienerSystem, A; kwargs...)
    length(sys.f) == 1 || error("describing_function requires a single nonlinearity, got $(length(sys.f))")
    isempty(sys.A) && sys.D == [0 1; 1 0] || error("Computing the describing function of a dynamic nonlinearity or a generic nonlinear system is not supported.")
    f = sys.f[1]
    if f isa Hysteresis
        # The system applies only the static callable of the Hysteresis struct; the
        # analytical overload describes the composite relay-with-hysteresis element
        # and would give the wrong answer here.
        describing_function_numerical(f, A; kwargs...)
    else
        describing_function(f, A; kwargs...)
    end
end

"""
    describing_function_plot(L, f; A_range, kwargs...)
    describing_function_plot(sys::HammersteinWienerSystem; A_range, kwargs...)

Create a Nyquist plot of the linear system `L` with the curve ``-1/N(A)`` overlaid,
where ``N(A)`` is the describing function of the nonlinearity `f`.

If a [`HammersteinWienerSystem`](@ref) is passed, the linear part and nonlinearity
are extracted automatically. The system must be an *open-loop* interconnection with a
single nonlinearity and a SISO external part, e.g., `G*saturation(1.0)`; the loop is
closed internally with unity negative feedback around the external channel. Do not pass
an already closed-loop system such as `feedback(G*nl)` — the loop would then be closed
twice and the plotted loop transfer would be wrong.

Intersections of the Nyquist curve and the ``-1/N(A)`` curve indicate potential
limit cycles. The ends of the ``-1/N(A)`` curve are marked to indicate the direction
of increasing amplitude; with the plotly backend, hovering over the curve shows the
amplitude ``A``.

# Arguments
- `L`: A linear system (transfer function or state-space).
- `f`: A static nonlinearity (callable or struct like [`Saturation`](@ref ControlSystemsBase.saturation)).
- `A_range`: Range of amplitudes for the describing function (default `range(0.1, 10, length=200)`).
- `kwargs...`: Additional keyword arguments passed to `nyquistplot`.

# Example
```julia
using ControlSystems, Plots
s = tf("s")
G = 10 / (s^3 + 2s^2 + s + 1)
describing_function_plot(G, Saturation(1.0); A_range=0.01:0.01:100)
```
"""
function describing_function_plot(L, f; A_range=range(0.1, 10, length=200), N=1000, kwargs...)
    neg_inv = map(A_range) do a
        n = describing_function(f, a; N)
        iszero(n) ? complex(NaN) : -1 / n
    end
    fig = nyquistplot(L; kwargs...)
    RecipesBase.plot!(fig, real.(neg_inv), imag.(neg_inv), lab="-1/N(A)", lw=2.5, arrow=true, hover = A_range)
    # Mark the curve endpoints to indicate the direction of increasing amplitude;
    # the arrow and hover information are not supported by all plotting backends.
    i1 = findfirst(isfinite, neg_inv)
    i2 = findlast(isfinite, neg_inv)
    if i1 !== nothing && i1 != i2
        RecipesBase.plot!(fig, [real(neg_inv[i1])], [imag(neg_inv[i1])], seriestype=:scatter, markershape=:circle, lab="A = $(round(A_range[i1], sigdigits=3))")
        RecipesBase.plot!(fig, [real(neg_inv[i2])], [imag(neg_inv[i2])], seriestype=:scatter, markershape=:square, lab="A = $(round(A_range[i2], sigdigits=3))")
    end
    fig
end

function describing_function_plot(sys::HammersteinWienerSystem; kwargs...)
    length(sys.f) == 1 || error("describing_function_plot requires a single nonlinearity, got $(length(sys.f))")
    sys.P.nu1 == 1 && sys.P.ny1 == 1 || error("describing_function_plot requires a SISO external part, got $(sys.P.ny1) external outputs and $(sys.P.nu1) external inputs.")
    # Close the outer (upper) loop with unity negative feedback to get the
    # linear transfer from the nonlinearity output (Δu) back to its input (Δy).
    # This is L_eff = -L_cl where L_cl is the positive-feedback loop transfer.
    # The condition for limit cycles is L_eff(jω) = -1/N_f(A), where N_f is
    # the DF of the actual static nonlinearity f (not a composite LFT system).
    L = sminreal(-lft(sys.P.P, ss(-1.0), :u))
    f = sys.f[1]
    f isa Hysteresis && @warn "This method of computing the describing function for a hysteresis element is known to be inaccurate. Hysteresis is a dynamic nonlinearity, and when approximated by an LFT system the linear part does not attenuate high frequencies as assumed by the describing function theory. For more accurate limit cycle predictions, consider using the `Hysteresis` struct (capital `H`) and call `describing_function_plot(P, H)` rather than forming the composite system by composition with `hysteresis(...)` (lowercase `h`). Computing the describing function with an instance of the `Hysteresis` struct uses the analytical formula for the describing function of an ideal hysteresis element."
    # Wrap in a closure to bypass analytical DF dispatch and force
    # describing_function_numerical (e.g., Hysteresis has an analytical DF for
    # the *composite* hysteresis element, but here we need the DF of the raw
    # static function since the internal linear dynamics are already included in L).
    describing_function_plot(L, x -> f(x); kwargs...)
end

