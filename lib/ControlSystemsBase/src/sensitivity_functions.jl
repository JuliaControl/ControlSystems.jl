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
PS = G[1:P.ny,     P.nu+1:end]
CS = G[P.ny+1:end, 1:P.nu]
T  = G[P.ny+1:end, P.nu+1:end]
```

The gang of four can be plotted like so
```julia
Gcl = extended_gangoffour(G, C) # Form closed-loop system
bodeplot(Gcl, lab=["S" "CS" "PS" "T"], plotphase=false) |> display # Plot gang of four
```
Note, the last output of Gcl is the negative of the `CS` and `T` transfer functions from `gangoffour2`. To get a transfer matrix with the same sign as [`G_CS`](@ref) and [`comp_sensitivity`](@ref), call `extended_gangoffour(P, C, pos=false)`.
See [`glover_mcfarlane`](@ref) from RobustAndOptimalControl.jl for an extended example. See also [`ncfmargin`](@ref) and [`feedback_control`](@ref) from RobustAndOptimalControl.jl.
"""
function extended_gangoffour(P, C, pos=true)
    ny,nu = size(P)
    te = P.timeevol
    if pos
        S = feedback(ss(I(ny+nu), P.timeevol), [ss(0*I(ny), te) P; -C ss(0*I(nu), te)], pos_feedback=true)
        return S + cat(0*I(ny), -I(nu), dims=(1,2))
    else
        Gtop = [I(ny); C] * [I(ny) P]
        return feedback(Gtop, ss(I(nu)), U1=(1:nu).+ny, Y1=(1:nu).+ny, pos_feedback=false)
    end
end

