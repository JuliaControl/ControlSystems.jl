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

