sensdoc = """
```
         ▲
         │e₁
         │  ┌─────┐
d₁────+──┴──►  P  ├─────┬──►e₄
      │     └─────┘     │
      │                 │
      │     ┌─────┐    -│
 e₂◄──┴─────┤  C  ◄──┬──+───d₂
            └─────┘  │
                     │e₃
                     ▼
```
- [`input_sensitivity`](@ref) is the transfer function from d₁ to e₁,       (I + CP)⁻¹
- [`output_sensitivity`](@ref) is the transfer function from d₂ to e₃,      (I + PC)⁻¹
- [`input_comp_sensitivity`](@ref) is the transfer function from d₁ to e₂,  (I + CP)⁻¹CP
- [`output_comp_sensitivity`](@ref) is the transfer function from d₂ to e₄, (I + PC)⁻¹PC
- [`G_PS`](@ref) is the transfer function from d₁ to e₄,                    (1 + PC)⁻¹P
- [`G_CS`](@ref) is the transfer function from d₂ to e₂,                    (1 + CP)⁻¹C
"""

"""
See [`output_sensitivity`](@ref)

The output [sensitivity function](https://en.wikipedia.org/wiki/Sensitivity_(control_systems)) ``S_o = (I + PC)^{-1}`` is the transfer function from a reference input to control error, while the input sensitivity function ``S_i = (I + CP)^{-1}`` is the transfer function from a disturbance at the plant input to the total plant input. For SISO systems, input and output sensitivity functions are equal. In general, we want to minimize the sensitivity function to improve robustness and performance, but practical constraints always cause the sensitivity function to tend to 1 for high frequencies. A robust design minimizes the peak of the sensitivity function, ``M_S``. The peak magnitude of ``S`` is the inverse of the distance between the open-loop Nyquist curve and the critical point -1. Upper bounding the sensitivity peak ``M_S`` gives lower-bounds on phase and gain margins according to
```math
ϕ_m ≥ 2\\text{sin}^{-1}(\\frac{1}{2M_S}), g_m ≥ \\frac{M_S}{M_S-1}
```
Generally, bounding ``M_S`` is a better objective than looking at gain and phase margins due to the possibility of combined gain and pahse variations, which may lead to poor robustness despite large gain and pahse margins.

$sensdoc
"""
function sensitivity(args...)# Sensitivity function
    return output_sensitivity(args...)
end

"""
See [`output_comp_sensitivity`](@ref)
$sensdoc
"""
function comp_sensitivity(args...) # Complementary sensitivity function
    return output_comp_sensitivity(args...)
end

"""
    G_PS(P, C)

The closed-loop transfer function from load disturbance to plant output.
Technically, the transfer function is given by `(1 + PC)⁻¹P` so `SP` would be a better, but nonstandard name.
$sensdoc
"""
G_PS(P, C) = feedback(P,C)#output_sensitivity(P, C)*P

"""
    G_CS(P, C)

The closed-loop transfer function from (-) measurement noise or (+) reference to control signal.
Technically, the transfer function is given by `(1 + CP)⁻¹C` so `SC` would be a better, but nonstandard name.
$sensdoc
"""
G_CS(P, C) = feedback(C,P)#input_sensitivity(P, C)*C

"""
    input_sensitivity(P, C)

Transfer function from load disturbance to total plant input.
- "Input" signifies that the transfer function is from the input of the plant.
$sensdoc
"""
function input_sensitivity(P,C)
    CP = C*P
    feedback(ss(Matrix{numeric_type(CP)}(I, ninputs(CP), ninputs(CP)), P.timeevol), CP)
end

"""
    output_sensitivity(P, C)

Transfer function from measurement noise / reference to control error.
- "output" signifies that the transfer function is from the output of the plant.
$sensdoc
"""
function output_sensitivity(P,C)
    PC = P*C
    S = feedback(ss(Matrix{numeric_type(PC)}(I, ninputs(PC), ninputs(PC)), P.timeevol), PC)
    S
end

"""
    input_comp_sensitivity(P,C)

Transfer function from load disturbance to control signal.
- "Input" signifies that the transfer function is from the input of the plant.
- "Complimentary" signifies that the transfer function is to an output (in this case controller output)
$sensdoc
"""
function input_comp_sensitivity(P,C)
    feedback(C * P)
end

"""
    output_comp_sensitivity(P,C)

Transfer function from measurement noise / reference to plant output.
- "output" signifies that the transfer function is from the output of the plant.
- "Complimentary" signifies that the transfer function is to an output (in this case plant output)
$sensdoc
"""
function output_comp_sensitivity(P,C)
    feedback(P * C)
end

"""
    extended_gangoffour(P, C, pos=true)

Returns a single statespace system that maps 
- `w1` reference or measurement noise
- `w2` load disturbance
to
- `z1` control error
- `z2` control input
```
      z1          z2
      ▲  ┌─────┐  ▲      ┌─────┐
      │  │     │  │      │     │
w1──+─┴─►│  C  ├──┴───+─►│  P  ├─┐
    │    │     │      │  │     │ │
    │    └─────┘      │  └─────┘ │
    │                 w2         │
    └────────────────────────────┘
```

The returned system has the transfer-function matrix
```math
\\begin{bmatrix}
I \\\\ C
\\end{bmatrix} (I + PC)^{-1} \\begin{bmatrix}
I & P
\\end{bmatrix}
```
or in code
```julia
# For SISO P
S  = G[1, 1]
PS = G[1, 2]
CS = G[2, 1]
T  = G[2, 2]

# For MIMO P
S  = G[1:P.ny,     1:P.nu]
PS = G[1:P.ny,     P.ny+1:end]
CS = G[P.ny+1:end, 1:P.ny]
T  = G[P.ny+1:end, P.ny+1:end] # Input complimentary sensitivity function
```

The gang of four can be plotted like so
```julia
Gcl = extended_gangoffour(G, C) # Form closed-loop system
bodeplot(Gcl, lab=["S" "CS" "PS" "T"], plotphase=false) |> display # Plot gang of four
```
Note, the last input of Gcl is the negative of the `PS` and `T` transfer functions from `gangoffour2`. To get a transfer matrix with the same sign as [`G_PS`](@ref) and [`input_comp_sensitivity`](@ref), call `extended_gangoffour(P, C, pos=false)`.
See [`glover_mcfarlane`](@ref) from RobustAndOptimalControl.jl for an extended example. See also [`ncfmargin`](@ref) and [`feedback_control`](@ref) from RobustAndOptimalControl.jl.
"""
function extended_gangoffour(P, C, pos=true)
    ny,nu = size(P)
    te = timeevol(P)
    if pos
        S = feedback(ss(I(ny+nu), te), [ss(0*I(ny), te) -P; C ss(0*I(nu), te)], pos_feedback=true)
        return S + cat(0*I(ny), -I(nu), dims=(1,2))
    else
        Gtop = [I(ny); C] * [I(ny) P]
        return feedback(Gtop, ss(I(nu)), U1=(1:nu).+ny, Y1=(1:nu).+ny, pos_feedback=false)
    end
end

