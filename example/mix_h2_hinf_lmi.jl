using LinearAlgebra
using JuMP
using Hypatia # You can choose other SDP solvers like Mosek, SCS, etc.

# Define the ExtendedStateSpace type with appropriate partitioning
# nx: number of states
# nu: number of control inputs
# ny: number of measured outputs
# nw_inf: number of H-infinity disturbance inputs
# nz_inf: number of H-infinity regulated outputs
# nw_2: number of H2 disturbance inputs
# nz_2: number of H2 regulated outputs
struct ExtendedStateSpace
    A::Matrix{Float64}
    B_w_inf::Matrix{Float64} # Bw_inf
    B_w_2::Matrix{Float64}   # Bw_2
    B_u::Matrix{Float64}     # Bu
    C_z_inf::Matrix{Float64} # Cz_inf
    D_z_inf_w_inf::Matrix{Float64} # Dz_inf_w_inf
    D_z_inf_u::Matrix{Float64}     # Dz_inf_u
    C_z_2::Matrix{Float64}   # Cz_2
    D_z_2_w_inf::Matrix{Float64}   # Dz_2_w_inf (often zero)
    D_z_2_w_2::Matrix{Float64}     # Dz_2_w_2 (must be zero for finite H2)
    D_z_2_u::Matrix{Float64}     # Dz_2_u
    C_y::Matrix{Float64}     # Cy (for state-feedback, this is identity for x, 0 for other, or unused)
    D_y_w_inf::Matrix{Float64} # Dyw_inf
    D_y_w_2::Matrix{Float64}   # Dyw_2
    D_y_u::Matrix{Float64}     # Dyu (often zero for state-feedback)

    nx::Int
    nu::Int
    nw_inf::Int
    nw_2::Int
    nz_inf::Int
    nz_2::Int
    ny::Int

    function ExtendedStateSpace(
        A, B_w_inf, B_w_2, B_u,
        C_z_inf, D_z_inf_w_inf, D_z_inf_u,
        C_z_2, D_z_2_w_inf, D_z_2_w_2, D_z_2_u,
        C_y, D_y_w_inf, D_y_w_2, D_y_u
    )
        nx = size(A, 1)
        nu = size(B_u, 2)
        nw_inf = size(B_w_inf, 2)
        nw_2 = size(B_w_2, 2)
        nz_inf = size(C_z_inf, 1)
        nz_2 = size(C_z_2, 1)
        ny = size(C_y, 1)

        # Basic dimension checks (can be made more exhaustive)
        @assert size(A) == (nx, nx)
        @assert size(B_w_inf) == (nx, nw_inf)
        @assert size(B_w_2) == (nx, nw_2)
        @assert size(B_u) == (nx, nu)

        @assert size(C_z_inf) == (nz_inf, nx)
        @assert size(D_z_inf_w_inf) == (nz_inf, nw_inf)
        @assert size(D_z_inf_u) == (nz_inf, nu)

        @assert size(C_z_2) == (nz_2, nx)
        @assert size(D_z_2_w_inf) == (nz_2, nw_inf)
        @assert size(D_z_2_w_2) == (nz_2, nw_2)
        @assert size(D_z_2_u) == (nz_2, nu)

        @assert size(C_y) == (ny, nx)
        @assert size(D_y_w_inf) == (ny, nw_inf)
        @assert size(D_y_w_2) == (ny, nw_2)
        @assert size(D_y_u) == (ny, nu)

        new(A, B_w_inf, B_w_2, B_u,
            C_z_inf, D_z_inf_w_inf, D_z_inf_u,
            C_z_2, D_z_2_w_inf, D_z_2_w_2, D_z_2_u,
            C_y, D_y_w_inf, D_y_w_2, D_y_u,
            nx, nu, nw_inf, nw_2, nz_inf, nz_2, ny
        )
    end
end

function robust_mixed_h2_hinf(
    systems::AbstractVector{<:ExtendedStateSpace};
    opt = Hypatia.Optimizer,
    gamma_inf_constraint::Float64, # The H-infinity bound to enforce
    verbose = true,
    silent_solver = true,
    ϵ = 1e-6, # Small positive constant for strict inequalities
)
    # Get dimensions from the first system (assume all systems have same dimensions)
    sys_nominal = systems[1]
    (; nx, nu, nw_inf, nw_2, nz_inf, nz_2) = sys_nominal

    model = JuMP.Model(opt)
    JuMP.set_optimizer_attribute(model, JuMP.MOI.Silent(), silent_solver)

    # LMI variables
    # X = P_inv (P is the common Lyapunov matrix)
    JuMP.@variable(model, X[1:nx, 1:nx], PSD)
    # Y = K * X (K is the state-feedback gain)
    JuMP.@variable(model, Y[1:nu, 1:nx])
    # Q_H2 for the H2 performance upper bound (its trace is minimized)
    JuMP.@variable(model, Q_H2[1:nz_2, 1:nz_2], PSD)

    # Objective: Minimize the trace of Q_H2 (upper bound on H2 norm squared)
    JuMP.@objective(model, Min, tr(Q_H2))

    # Add constraints for each plant model in the set (polytopic vertices)
    for (idx, sys) in enumerate(systems)
        # Check for D_z_2_w_2 == 0 for H2 finite norm
        if !iszero(sys.D_z_2_w_2)
            error("D_z_2_w_2 must be zero for finite H2 norm in this formulation for system $idx.")
        end

        # H-infinity Constraint LMI (for each vertex)
        # [ A_i X + X A_i^T + B_ui Y + Y^T B_ui^T   B_wi   X C_zi^T + Y^T D_zui^T ]
        # [ B_wi^T                                  -gamma_inf I  D_zwi^T            ]
        # [ C_zi X + D_zui Y                      D_zwi      -gamma_inf I         ]

        # Note: B_wi here refers to B_w_inf, and so on for C_zi, D_zwi, D_zui
        # To avoid confusion, I'll use full variable names from the struct
        Hinf_LMI_block11 = sys.A * X + X * sys.A' + sys.B_u * Y + Y' * sys.B_u'
        Hinf_LMI_block12 = sys.B_w_inf
        Hinf_LMI_block13 = X * sys.C_z_inf' + Y' * sys.D_z_inf_u'

        Hinf_LMI_block21 = sys.B_w_inf'
        Hinf_LMI_block22 = -gamma_inf_constraint * I(nw_inf)
        Hinf_LMI_block23 = sys.D_z_inf_w_inf' # This term is often 0 or simplified if D_z_inf_w_inf is 0

        Hinf_LMI_block31 = sys.C_z_inf * X + sys.D_z_inf_u * Y
        Hinf_LMI_block32 = sys.D_z_inf_w_inf
        Hinf_LMI_block33 = -gamma_inf_constraint * I(nz_inf)

        Hinf_LMI = [
            Hinf_LMI_block11 Hinf_LMI_block12 Hinf_LMI_block13;
            Hinf_LMI_block21 Hinf_LMI_block22 Hinf_LMI_block23;
            Hinf_LMI_block31 Hinf_LMI_block32 Hinf_LMI_block33;
        ]
        JuMP.@constraint(model, Hinf_LMI <= -ϵ * I, PSDCone())

        # H2 Performance LMI (for each vertex)
        # [ A_i X + X A_i^T + B_ui Y + Y^T B_ui^T   B_w2i   X C_z2i^T + Y^T D_z2ui^T ]
        # [ B_w2i^T                                  -I      0                        ]
        # [ C_z2i X + D_z2ui Y                      0       -Q_H2                    ]

        H2_LMI_block11 = sys.A * X + X * sys.A' + sys.B_u * Y + Y' * sys.B_u'
        H2_LMI_block12 = sys.B_w_2
        H2_LMI_block13 = X * sys.C_z_2' + Y' * sys.D_z_2_u'

        H2_LMI_block21 = sys.B_w_2'
        H2_LMI_block22 = -I(nw_2)
        H2_LMI_block23 = zeros(nw_2, nz_2) # Assuming D_z_2_w_2 = 0

        H2_LMI_block31 = sys.C_z_2 * X + sys.D_z_2_u * Y
        H2_LMI_block32 = zeros(nz_2, nw_2) # Assuming D_z_2_w_2 = 0
        H2_LMI_block33 = -Q_H2

        H2_LMI = [
            H2_LMI_block11 H2_LMI_block12 H2_LMI_block13;
            H2_LMI_block21 H2_LMI_block22 H2_LMI_block23;
            H2_LMI_block31 H2_LMI_block32 H2_LMI_block33;
        ]
        JuMP.@constraint(model, H2_LMI <= -ϵ * I, PSDCone())
    end

    # Explicitly add X positive definite constraint if not implicitly covered by the large matrices
    JuMP.@constraint(model, X >= ϵ * I, PSDCone())


    JuMP.optimize!(model)

    if verbose
        @info "JuMP Termination Status: $(JuMP.termination_status(model))"
        if JuMP.has_values(model)
            @info "JuMP Primal Status: $(JuMP.primal_status(model))"
            @info "Guaranteed H2 Cost (tr(Q_H2)): $(JuMP.objective_value(model))"
        else
            @warn "No solution found."
        end
    end

    if JuMP.has_values(model) && JuMP.termination_status(model) == JuMP.MOI.OPTIMAL
        X_val = JuMP.value.(X)
        Y_val = JuMP.value.(Y)
        Q_H2_val = JuMP.value.(Q_H2)

        # Calculate K = Y * X_inv
        K = Y_val * inv(X_val)
        guaranteed_h2_cost = tr(Q_H2_val)

        return (K=K, X=X_val, Q_H2=Q_H2_val, guaranteed_h2_cost=guaranteed_h2_cost)
    else
        return (K=nothing, X=nothing, Q_H2=nothing, guaranteed_h2_cost=nothing)
    end
end


## =============================================================================
# Example Usage:
# Define dimensions
nx = 2  # States
nu = 1  # Control inputs
nw_inf = 1 # H-infinity disturbance inputs
nw_2 = 1   # H2 disturbance inputs
nz_inf = 2 # H-infinity regulated outputs (e.g., weighted x[1] and u)
nz_2 = 2   # H2 regulated outputs (e.g., weighted x and u)
ny = nx    # Measured outputs (assuming state-feedback, so all states measured)

# Define nominal system matrices (Vertex 1)
A1 = [-1.0 0.5; 0.2 -2.0]
Bu1 = [1.0; 0.5;;]
Bw_inf1 = [0.1; 0.05;;]
Bw_21 = [0.0; 0.1;;] # For H2, usually uncorrelated disturbances

# Performance outputs for H-infinity
# Let z_inf = [x_1; 0.1*u]
Cz_inf1 = [1.0 0.0; 0.0 0.0]
Dz_inf_u1 = [0.0; 0.1;;]
Dz_inf_w_inf1 = zeros(nz_inf, nw_inf) # Often zero for H-inf performance

# Performance outputs for H2
# Let z_2 = [0.5*x; 0.1*u]
Cz_21 = [0.5 0.0; 0.0 0.5]
Dz_2_u1 = [0.0; 0.1;;]
Dz_2_w_inf1 = zeros(nz_2, nw_inf) # H2 part does not depend on H-inf disturbances directly
Dz_2_w_21 = zeros(nz_2, nw_2) # Must be zero for finite H2 norm

# Measured output for state-feedback (all states measured)
Cy1 = I(nx)
Dyw_inf1 = zeros(ny, nw_inf)
Dyw_21 = zeros(ny, nw_2)
Dyu1 = zeros(ny, nu)

sys1 = ExtendedStateSpace(
    A1, Bw_inf1, Bw_21, Bu1,
    Cz_inf1, Dz_inf_w_inf1, Dz_inf_u1,
    Cz_21, Dz_2_w_inf1, Dz_2_w_21, Dz_2_u1,
    Cy1, Dyw_inf1, Dyw_21, Dyu1
)

# Define uncertain system matrices (Vertex 2)
# Slightly different A matrix to represent uncertainty
A2 = [-1.2 0.4; 0.1 -2.2]
Bu2 = [1.1; 0.4;;]
Bw_inf2 = [0.1; 0.05;;] # Same disturbances for simplicity
Bw_22 = [0.0; 0.1;;]

# Same performance output definitions for the second vertex
Cz_inf2 = Cz_inf1
Dz_inf_u2 = Dz_inf_u1
Dz_inf_w_inf2 = Dz_inf_w_inf1

Cz_22 = Cz_21
Dz_2_u2 = Dz_2_u1
Dz_2_w_inf2 = Dz_2_w_inf1
Dz_2_w_22 = Dz_2_w_21

Cy2 = Cy1
Dyw_inf2 = Dyw_inf1
Dyw_22 = Dyw_21
Dyu2 = Dyu1

sys2 = ExtendedStateSpace(
    A2, Bw_inf2, Bw_22, Bu2,
    Cz_inf2, Dz_inf_w_inf2, Dz_inf_u2,
    Cz_22, Dz_2_w_inf2, Dz_2_w_22, Dz_2_u2,
    Cy2, Dyw_inf2, Dyw_22, Dyu2
)

# Array of uncertain plant models (vertices of the polytope)
plant_models = [sys1, sys2]

# Desired H-infinity constraint
gamma_inf_desired = 0.5 # Try to achieve an H-inf norm less than 0.5

# Solve the mixed H2/H-infinity problem
result = robust_mixed_h2_hinf(
    plant_models,
    gamma_inf_constraint = gamma_inf_desired,
    verbose = true,
    silent_solver = false # Set to false to see solver output
)

if result.K !== nothing
    println("\nController Gain K:")
    display(result.K)
    println("\nCommon Lyapunov Matrix Inverse (X):")
    display(result.X)
    println("\nGuaranteed H2 Cost (trace(Q_H2)): ", result.guaranteed_h2_cost)
else
    println("\nFailed to find a feasible solution for the given constraints.")
end