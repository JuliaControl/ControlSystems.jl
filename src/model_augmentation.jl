"""
    add_disturbance(sys::AbstractStateSpace{Continuous}, Ad::AbstractMatrix, Cd::AbstractMatrix)

See CCS pp. 144

# Arguments:
- `sys`: System to augment
- `Ad`: The dynamics of the disturbance
- `Cd`: How the disturbance states affect the states of `sys`. This matrix as the shape (sys.nx, size(Ad, 1))

See also `add_low_frequency_disturbance, add_resonant_disturbance`
"""
function add_disturbance(sys::AbstractStateSpace{Continuous}, Ad::AbstractMatrix, Cd::AbstractMatrix)
    A,B,C,D = ControlSystems.ssdata(sys)
    T = eltype(A)
    nx,nu,ny = sys.nx,sys.nu,sys.ny
    Ae = [A Cd; zeros(T, size(Ad, 1), nx) Ad]
    Be = [B; zeros(T, size(Ad, 1), nu)]
    Ce = [C zeros(T, ny, size(Ad, 1))]
    De = D
    ss(Ae,Be,Ce,De)
end

function add_measurement_disturbance(sys::AbstractStateSpace{Continuous}, Ad::AbstractMatrix, Cd::AbstractMatrix)
    A,B,C,D = ControlSystems.ssdata(sys)
    T = eltype(A)
    nx,nu,ny = sys.nx,sys.nu,sys.ny
    Ae = [A zeros(T, nx, size(Ad, 1)); zeros(T, size(Ad, 1), nx) Ad]
    Be = [B; zeros(T, size(Ad, 1), nu)]
    Ce = [C Cd]
    De = D
    ss(Ae,Be,Ce,De)
end

function add_low_frequency_disturbance(sys::AbstractStateSpace{Continuous}, Ai::Integer; ϵ=0)
    nx,nu,ny = sys.nx,sys.nu,sys.ny
    Cd = zeros(nx, 1)
    Cd[Ai] = 1
    add_disturbance(sys, fill(-ϵ, 1, 1), Cd)
end

function add_low_frequency_disturbance(sys::AbstractStateSpace{Continuous}; ϵ=0, measurement=false)
    nx,nu,ny = sys.nx,sys.nu,sys.ny
    if measurement
        Cd = I(nu)
        add_measurement_disturbance(sys, -ϵ*I(nu), Cd)
    else
        Cd = sys.B
        add_disturbance(sys, -ϵ*I(nu), Cd)
    end
end

function add_resonant_disturbance(sys::AbstractStateSpace{Continuous}, ω, ζ, Ai::Integer; measurement=false)
    nx,nu,ny = sys.nx,sys.nu,sys.ny
    if measurement
        Cd = zeros(ny, 2)
        Cd[Ai, 1] = 1
    else
        Cd = zeros(nx, 2)
        Cd[Ai, 1] = 1
    end
    Ad = [-ζ -ω; ω -ζ]
    measurement ? add_measurement_disturbance(sys, Ad, Cd) : add_disturbance(sys, Ad, Cd)
end
