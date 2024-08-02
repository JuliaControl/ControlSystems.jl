#####################################################################
##                      Data Type Declarations                     ##
#####################################################################

function state_space_validation(A,B,C,D)
    nx = size(A, 1)
    nu = size(B, 2)
    ny = size(C, 1)

    if size(A, 2) != nx && nx != 0
        error("A has dimensions $(size(A)), but must be square")
    elseif size(B, 1) != nx
        error("The number of rows of A ($(size(A,1))) and B ($(size(B,1))) are not equal")
    elseif size(C, 2) != nx
        error("The number of columns of A ($(size(A,2))) and C ($(size(C,2))) are not equal")
    elseif nu != size(D, 2)
        error("The number of columns of B ($(size(B,2))) and D ($(size(D,2))) are not equal")
    elseif ny != size(D, 1)
        error("The number of rows of C ($(size(C,1))) and D ($(size(D,1))) are not equal")
    end

    nx,nu,ny
end

abstract type AbstractStateSpace{TE<:TimeEvolution} <: LTISystem{TE} end

"""
    StateSpace{TE, T} <: AbstractStateSpace{TE}

An object representing a standard state space system.

See the function [`ss`](@ref) for a user facing constructor as well as the documentation page [creating systems](https://juliacontrol.github.io/ControlSystems.jl/stable/man/creating_systems/).

# Fields:
- `A::Matrix{T}`
- `B::Matrix{T}`
- `C::Matrix{T}`
- `D::Matrix{T}`
- `timeevol::TE`
"""
struct StateSpace{TE, T} <: AbstractStateSpace{TE}
    A::Matrix{T}
    B::Matrix{T}
    C::Matrix{T}
    D::Matrix{T}
    timeevol::TE
    function StateSpace{TE, T}(A::Matrix{T}, B::Matrix{T}, C::Matrix{T}, D::Matrix{T}, timeevol::TE) where {TE<:TimeEvolution, T}
        state_space_validation(A,B,C,D)
        new{TE, T}(A, B, C, D, timeevol)
    end
    function StateSpace(A::Matrix{T}, B::Matrix{T}, C::Matrix{T}, D::Matrix{T}, timeevol::TE) where {TE<:TimeEvolution, T}
        state_space_validation(A,B,C,D)
        new{TE, T}(A, B, C, D, timeevol)
    end
end
# Constructor for Discrete system
function StateSpace(A::Matrix{T}, B::Matrix{T}, C::Matrix{T}, D::Matrix{T}, Ts::Number) where {T}
    state_space_validation(A,B,C,D)
    StateSpace(A, B, C, D, Discrete(Ts))
end
# Constructor for Continuous system
function StateSpace(A::Matrix{T}, B::Matrix{T}, C::Matrix{T}, D::Matrix{T}) where {T}
    state_space_validation(A,B,C,D)
    StateSpace(A, B, C, D, Continuous())
end

""" If D=0 then create zero matrix of correct size and type, else, convert D to correct type"""
function fix_D_matrix(T::Type,B,C,D)
    if D == 0
        D = zeros(T, size(C,1), size(B,2))
    else
        D = to_matrix(T, D)
    end
    return D
end
""" Promote A,B,C,D to same types, fix D matrix (see `fix_D_matrix`)"""
function to_similar_matrices(A,B,C,D)
    T = promote_type(eltype(A),eltype(B),eltype(C),eltype(D))
    A = to_matrix(T, A)
    B = to_matrix(T, B)
    C = C isa UniformScaling ? to_matrix(T, Matrix(C(size(A, 1)))) : to_matrix(T, C)
    D = fix_D_matrix(T, B, C, D)
    return A, B, C, D, T
end

const AbstractNumOrArray = Union{Number, AbstractVecOrMat}

# Explicit construction with different types
function StateSpace{TE,T}(A::AbstractNumOrArray, B::AbstractNumOrArray, C::AbstractNumOrArray, D::AbstractNumOrArray, timeevol::TimeEvolution) where {TE, T}
    D = fix_D_matrix(T, B, C, D)
    return StateSpace{TE, T}(to_matrix(T, A), to_matrix(T, B), to_matrix(T, C), to_matrix(T, D), TE(timeevol))
end

# Explicit conversion
function StateSpace{TE,T}(sys::StateSpace) where {TE, T}
    StateSpace{TE,T}(sys.A,sys.B,sys.C,sys.D,sys.timeevol)
end

function StateSpace(A::AbstractNumOrArray, B::AbstractNumOrArray, C::Union{AbstractNumOrArray, UniformScaling}, D::AbstractNumOrArray, timeevol::TimeEvolution)
    A, B, C, D, T = to_similar_matrices(A,B,C,D)
    return StateSpace{typeof(timeevol),T}(A, B, C, D, timeevol)
end
# General Discrete constructor
StateSpace(A::AbstractNumOrArray, B::AbstractNumOrArray, C::Union{AbstractNumOrArray, UniformScaling}, D::AbstractNumOrArray, Ts::Number) =
    StateSpace(A, B, C, D, Discrete(Ts))
# General continuous constructor
StateSpace(A::AbstractNumOrArray, B::AbstractNumOrArray, C::Union{AbstractNumOrArray, UniformScaling}, D::AbstractNumOrArray) =
    StateSpace(A, B, C, D, Continuous())

# Function for creation of static gain
function StateSpace(D::AbstractArray{T}, timeevol::TimeEvolution) where {T<:Number}
    ny, nu = size(D, 1), size(D, 2)
    A = zeros(T, 0, 0)
    B = zeros(T, 0, nu)
    C = zeros(T, ny, 0)
    D = reshape(D, (ny,nu))
    return StateSpace(A, B, C, D, timeevol)
end
StateSpace(D::AbstractArray, Ts::Number) = StateSpace(D, Discrete(Ts))
StateSpace(D::AbstractArray) = StateSpace(D, Continuous())

StateSpace(d::Number, Ts::Number; kwargs...) = StateSpace([d], Discrete(Ts))
StateSpace(d::Number; kwargs...) = StateSpace([d], Continuous())


# StateSpace(sys) converts to StateSpace
StateSpace(sys::LTISystem; kwargs...) = convert(StateSpace, sys; kwargs...)

"""
    sys = ss(A, B, C, D)      # Continuous
    sys = ss(A, B, C, D, Ts)  # Discrete

Create a state-space model `sys::StateSpace{TE, T}`
with matrix element type `T` and TE is `Continuous` or `<:Discrete`.

This is a continuous-time model if `Ts` is omitted.
Otherwise, this is a discrete-time model with sampling period `Ts`.

`D` may be specified as `0` in which case a zero matrix of appropriate size is constructed automatically. 
`sys = ss(D [, Ts])` specifies a static gain matrix `D`.

To associate names with states, inputs and outputs, see [`named_ss`](https://juliacontrol.github.io/RobustAndOptimalControl.jl/dev/#Named-systems).
"""
ss(args...;kwargs...) = StateSpace(args...;kwargs...)


struct HeteroStateSpace{TE <: TimeEvolution, AT<:AbstractMatrix,BT<:AbstractVecOrMat,CT<:AbstractMatrix,DT<:AbstractMatrix} <: AbstractStateSpace{TE}
    A::AT
    B::BT
    C::CT
    D::DT
    timeevol::TE
    function HeteroStateSpace(A::AT, B::BT,
        C::CT, D::DT, timeevol::TE) where {TE<:TimeEvolution,AT<:AbstractMatrix,BT<:AbstractMatrix,CT<:AbstractMatrix,DT<:AbstractMatrix}
        state_space_validation(A,B,C,D)
        new{TE,AT,BT,CT,DT}(A, B, C, D, timeevol)
    end
    # Explicit constructor
    function HeteroStateSpace{TE,AT,BT,CT,DT}(A, B, C, D, timeevol) where {TE,AT,BT,CT,DT}
        state_space_validation(A,B,C,D)
        new{TE,AT,BT,CT,DT}(
            A isa AT ? A : AT(A),
            B isa BT ? B : BT(B),
            C isa CT ? C : CT(C),
            D isa DT ? D : DT(D),
            TE(timeevol)
        )
    end
    # Base constructor
    function HeteroStateSpace(A::AbstractNumOrArray, B::AbstractNumOrArray, C::Union{AbstractNumOrArray, UniformScaling}, D::AbstractNumOrArray, timeevol::TimeEvolution)
        A = to_abstract_matrix(A)
        B = to_abstract_matrix(B)
        C = C isa UniformScaling ? C(size(A, 1)) : to_abstract_matrix(C)
        if D == 0
            D = fill(zero(eltype(C)), size(C,1), size(B,2))
        else
            D = to_abstract_matrix(D)
        end
        state_space_validation(A,B,C,D)
        return new{typeof(timeevol),typeof(A),typeof(B),typeof(C),typeof(D)}(A, B, C, D, timeevol)
    end
end

function HeteroStateSpace(A::AbstractNumOrArray, B::AbstractNumOrArray, C::Union{AbstractNumOrArray, UniformScaling}, D::AbstractNumOrArray, Ts::Number)
    HeteroStateSpace(A, B, C, D, Discrete(Ts))
end
function HeteroStateSpace(A::AbstractNumOrArray, B::AbstractNumOrArray, C::Union{AbstractNumOrArray, UniformScaling}, D::AbstractNumOrArray)
    HeteroStateSpace(A, B, C, D, Continuous())
end
HeteroStateSpace(s::AbstractStateSpace) = HeteroStateSpace(s.A,s.B,s.C,s.D,s.timeevol)



# Function for creation of static gain
function HeteroStateSpace(D::AbstractArray{T}, timeevol::TimeEvolution) where {T<:Number}
    ny, nu = size(D, 1), size(D, 2)
    A = zeros(T, 0, 0)
    B = zeros(T, 0, nu)
    C = zeros(T, ny, 0)

    return HeteroStateSpace(A, B, C, D, timeevol)
end

HeteroStateSpace(D::AbstractArray{T}, Ts::Number) where {T<:Number} = HeteroStateSpace(D, Discrete(Ts))
HeteroStateSpace(D::AbstractArray{T}) where {T<:Number} = HeteroStateSpace(D, Continuous())

HeteroStateSpace(d::Number, Ts::Number; kwargs...) = HeteroStateSpace([d], Discrete(Ts))
HeteroStateSpace(d::Number; kwargs...) = HeteroStateSpace([d], Continuous())

# HeteroStateSpace(sys) converts to HeteroStateSpace
HeteroStateSpace(sys::LTISystem) = convert(HeteroStateSpace, sys)

# Getter functions ###################################################

"""
    A, B, C, D = ssdata(sys)

A destructor that outputs the statespace matrices.
"""
ssdata(sys) = sys.A, sys.B, sys.C, sys.D

# Functions for number of inputs, outputs and states
ninputs(sys::AbstractStateSpace) = size(sys.D, 2)
noutputs(sys::AbstractStateSpace) = size(sys.D, 1)
nstates(sys::AbstractStateSpace) = size(sys.A, 1)

isproper(sys::AbstractStateSpace) = iszero(sys.D)

function isstable(sys::StateSpace{Continuous, <:ForwardDiff.Dual})
    # Drop duals for this check since it's not differentiable anyway
    all(real.(eigvals(ForwardDiff.value.(sys.A))) .< 0)
end
function isstable(sys::StateSpace{<:Discrete, <:ForwardDiff.Dual})
    all(abs.(ForwardDiff.value.(sys.A)) .< 1)
end


#####################################################################
##                         Math Operators                          ##
#####################################################################

## EQUALITY ##
function ==(sys1::ST1, sys2::ST2) where {ST1<:AbstractStateSpace,ST2<:AbstractStateSpace}
    fieldnames(ST1) == fieldnames(ST2) || (return false)
    return all(getfield(sys1, f) == getfield(sys2, f) for f in fieldnames(ST1))
end

## Approximate ##
function isapprox(sys1::ST1, sys2::ST2; kwargs...) where {ST1<:AbstractStateSpace,ST2<:AbstractStateSpace}
    fieldnames(ST1) == fieldnames(ST2) || (return false)
    return all(isapprox(getfield(sys1, f), getfield(sys2, f); kwargs...) for f in fieldnames(ST1))
end

## ADDITION ##
Base.zero(sys::AbstractStateSpace) = basetype(sys)(zero(sys.D), sys.timeevol)
Base.zero(::Type{StateSpace{Continuous, F}}) where {F} = ss([zero(F)], Continuous())
Base.zero(::Type{StateSpace{D, F}}) where {D<:Discrete, F} = ss([zero(F)], undef_sampletime(D))

function +(s1::ST, s2::ST) where {ST <: AbstractStateSpace}
    #Ensure systems have same dimensions
    if size(s1) != size(s2)
        error("Systems have different shapes.")
    end
    timeevol = common_timeevol(s1,s2)
    T = promote_type(numeric_type(s1), numeric_type(s2))

    A = [s1.A                   zeros(T, nstates(s1), nstates(s2));
         zeros(T, nstates(s2), nstates(s1))        s2.A]
    B = [s1.B ; s2.B]
    C = [s1.C s2.C;]
    D = [s1.D + s2.D;]

    return ST(A, B, C, D, timeevol)
end


function +(sys::ST, n::Number) where ST <: AbstractStateSpace
    issiso(sys) || throw(DimensionMismatch("Numbers and systems can only be added for SISO systems."))
    basetype(ST)(sys.A, sys.B, sys.C, sys.D .+ n, sys.timeevol)
end
+(n::Number, sys::ST) where ST <: AbstractStateSpace = +(sys, n)

## SUBTRACTION ##
-(sys1::AbstractStateSpace, sys2::AbstractStateSpace) = +(sys1, -sys2)
-(sys::AbstractStateSpace, n::Number) = +(sys, -n)
-(n::Number, sys::AbstractStateSpace) = +(-sys, n)

## NEGATION ##
-(sys::ST) where ST <: AbstractStateSpace = basetype(ST)(sys.A, sys.B, -sys.C, -sys.D, sys.timeevol)

## MULTIPLICATION ##
function *(sys1::ST, sys2::ST) where {ST <: AbstractStateSpace}
    #Check dimension alignment
    #Note: sys1*sys2 = y <- sys1 <- sys2 <- u
    if (sys1.nu != sys2.ny) && (sys1.nu == 1 || sys2.ny == 1)
        throw(DimensionMismatch("sys1*sys2: sys1 must have same number of inputs as sys2 has outputs. If you want to broadcast a scalar system to a diagonal system, use broadcasted multiplication sys1 .* sys2"))
    end
    sys1.nu == sys2.ny || throw(DimensionMismatch("sys1*sys2: sys1 must have same number of inputs as sys2 has outputs"))
    timeevol = common_timeevol(sys1,sys2)
    T = promote_type(numeric_type(sys1), numeric_type(sys2))

    A = [sys1.A    sys1.B*sys2.C;
         zeros(T, sys2.nx, sys1.nx)  sys2.A]
    B = [sys1.B*sys2.D ; sys2.B]
    C = [sys1.C   sys1.D*sys2.C;]
    D = [sys1.D*sys2.D;]
    return basetype(ST)(A, B, C, D, timeevol)
end

function Base.Broadcast.broadcasted(::typeof(*), sys1::AbstractStateSpace, sys2::AbstractStateSpace)
    issiso(sys1) || issiso(sys2) || throw(DimensionMismatch("Only SISO statespace systems can be broadcasted"))
    if issiso(sys1) && !issiso(sys2) # Check !issiso(sys2) to avoid calling fill if both are siso
        sys1 = append(sys1 for i in 1:sys2.ny)
    elseif issiso(sys2)
        sys2 = append(sys2 for i in 1:sys1.nu)
    end
    return sys1 * sys2
end

function Base.Broadcast.broadcasted(::typeof(*), sys1::ST, M::AbstractArray{<:Number}) where {ST <: AbstractStateSpace}
    LinearAlgebra.isdiag(M) || throw(DimensionMismatch("Broadcasting multiplication of an LTI system with an array is only supported for diagonal arrays. If you want the system to behave like a scalar and multiply each element of the array, wrap the system in a `Ref` to indicate this, i.e., `Ref(sys) .* array`. See also function `array2mimo`."))
    sys1 .* ss(M, sys1.timeevol) # If diagonal, broadcast by replicating input channels
end

function Base.Broadcast.broadcasted(::typeof(*), M::AbstractArray{<:Number}, sys1::ST) where {ST <: AbstractStateSpace}
    LinearAlgebra.isdiag(M) || throw(DimensionMismatch("Broadcasting multiplication of an LTI system with an array is only supported for diagonal arrays. If you want the system to behave like a scalar and multiply each element of the array, wrap the system in a `Ref` to indicate this, i.e., `array .* Ref(sys)`. See also function `array2mimo`."))
    ss(M, sys1.timeevol) .* sys1 # If diagonal, broadcast by replicating output channels
end

function Base.Broadcast.broadcasted(::typeof(*), sys1::Base.RefValue{ST}, M::AbstractArray{<:Number}) where {ST <: AbstractStateSpace}
    sys1 = sys1[]
    issiso(sys1) || throw(DimensionMismatch("Only SISO statespace systems can be broadcasted"))
    T = promote_type(numeric_type(sys1), eltype(M))
    A,B,C,D = ssdata(sys1)
    nx = sys1.nx
    ny,nu = size(M,1), size(M,2)
    Ae = cat(fill(A, ny*nu)..., dims=(1,2))
    Be::Matrix{T} = ControlSystemsBase.blockdiag(B for i in 1:nu)
    Be = repeat(Be, ny, 1)

    Ce::Matrix{T} = repeat(C, 1, nu)
    Ce = ControlSystemsBase.blockdiag(Ce for i in 1:ny)
    Ce = Ce * kron(Diagonal(vec(M')), I(nx))

    De = D .* M
    sminreal(basetype(ST)(Ae, Be, Ce, De, sys1.timeevol))
end

function Base.Broadcast.broadcasted(::typeof(*), M::AbstractArray{<:Number}, sys1::Base.RefValue{ST}) where {ST <: AbstractStateSpace}
    sys1 = sys1[]
    issiso(sys1) || throw(DimensionMismatch("Only SISO statespace systems can be broadcasted"))
    T = promote_type(numeric_type(sys1), eltype(M))
    A,B,C,D = ssdata(sys1)
    nx = sys1.nx
    ny,nu = size(M,1), size(M,2)
    Ae = cat(fill(A, ny*nu)..., dims=(1,2))
    Be::Matrix{T} = ControlSystemsBase.blockdiag(convert(Matrix{T}, B) for i in 1:nu)
    Be = repeat(Be, ny, 1)
    Be = kron(Diagonal(vec(M')), I(nx)) * Be

    Ce::Matrix{T} = repeat(convert(Matrix{T}, C), 1, nu)
    Ce = ControlSystemsBase.blockdiag(Ce for i in 1:ny)

    De = D .* M
    sminreal(basetype(ST)(Ae, Be, Ce, De, sys1.timeevol))
end

function *(sys1::ST, D::Diagonal) where {ST <: AbstractStateSpace}
    if issiso(sys1) # This is a special case that falls back on broadcasting
        return sys1 .* D
    else # This is the standard implementation but must be handled here since we special case diagonal matrices for the case above
        sys1 * ss(D, sys1.timeevol)
    end
end

*(sys::ST, n::Number) where ST <: AbstractStateSpace = basetype(ST)(sys.A, sys.B*n, sys.C, sys.D*n, sys.timeevol)
*(n::Number, sys::ST) where ST <: AbstractStateSpace = basetype(ST)(sys.A, sys.B, sys.C*n, sys.D*n, sys.timeevol)


## DIVISION ##


"""
    /(sys1::AbstractStateSpace{TE}, sys2::AbstractStateSpace{TE}; atol::Real = 0, atol1::Real = atol, atol2::Real = atol, rtol::Real = max(size(sys1.A, 1), size(sys2.A, 1)) * eps(real(float(one(numeric_type(sys1))))) * iszero(min(atol1, atol2)))

Compute `sys1 / sys2 = sys1 * inv(sys2)` in a way that tries to handle situations in which the inverse `sys2` is non-proper, but the resulting system `sys1 / sys2` is proper.

See `ControlSystemsBase.MatrixPencils.isregular` for keyword arguments `atol`, `atol1`, `atol2`, and `rtol`.
"""
function Base.:(/)(sys1::AbstractStateSpace{TE}, sys2::AbstractStateSpace{TE}; 
    atol::Real = 0, atol1::Real = atol, atol2::Real = atol, 
    rtol::Real = max(size(sys1.A,1),size(sys2.A,1))*eps(real(float(one(numeric_type(sys1)))))*iszero(min(atol1,atol2))) where {TE<:ControlSystemsBase.TimeEvolution}
    T1 = float(numeric_type(sys1))
    T2 = float(numeric_type(sys2))
    T = promote_type(T1,T2)
    timeevol = common_timeevol(sys1, sys2)
    ny2, nu2 = sys2.ny, sys2.nu
    nu2 == ny2  || error("The system sys2 must be square")
    ny1, nu1 = sys1.ny, sys1.nu
    nu1 == nu2  || error("The systems sys1 and sys2 must have the same number of inputs")
    nx1 = sys1.nx
    nx2 = sys2.nx
    if nx2 > 0
        A, B, C, D = ssdata([sys2; sys1])
        Ai = [A B; C[1:ny2,:] D[1:ny2,:]]
        Ei = [I zeros(T,nx1+nx2,ny2); zeros(T,ny2,nx1+nx2+ny2)] |> Matrix # TODO: rm call to Matrix when type piracy in https://github.com/JuliaLinearAlgebra/LinearMaps.jl/issues/219 is fixed
        MatrixPencils.isregular(Ai, Ei; atol1, atol2, rtol) || 
            error("The system sys2 is not invertible")
        Ci = [C[ny2+1:ny1+ny2,:] D[ny2+1:ny1+ny2,:]]
        Bi = [zeros(T,nx1+nx2,nu1); -I] |> Matrix # TODO: rm call to Matrix when type piracy in https://github.com/JuliaLinearAlgebra/LinearMaps.jl/issues/219 is fixed
        Di = zeros(T,ny1,nu1)
        Ai, Ei, Bi, Ci, Di = MatrixPencils.lsminreal(Ai, Ei, Bi, Ci, Di; fast = true, atol1 = 0, atol2, rtol, contr = true, obs = true, noseig = true)
        if Ei != I
            luE = lu!(Ei, check=false)
            issuccess(luE) || throw(ArgumentError("The system sys2 is not invertible"))
            Ai = luE\Ai
            Bi = luE\Bi
        end
    else
        D2 = T.(sys2.D)
        LUD = lu(D2)
        (norm(D2,Inf) <= atol1 || rcond(LUD.U) <= 10*nu1*eps(real(float(one(T))))) && 
                error("The system sys2 is not invertible")
        Ai, Bi, Ci, Di = ssdata(sys1)
        rdiv!(Bi,LUD); rdiv!(Di,LUD)
    end

    return StateSpace{TE, T}(Ai, Bi, Ci, Di, timeevol) 
end

function Base.:(\)(sys1::AbstractStateSpace{TE}, sys2::AbstractStateSpace{TE}; 
    atol::Real = 0, atol1::Real = atol, atol2::Real = atol, 
    rtol::Real = max(size(sys1.A,1),size(sys2.A,1))*eps(real(float(one(numeric_type(sys1)))))*iszero(min(atol1,atol2))) where {TE<:ControlSystemsBase.TimeEvolution}
    T1 = float(numeric_type(sys1))
    T2 = float(numeric_type(sys2))
    T = promote_type(T1,T2)
    timeevol = common_timeevol(sys1, sys2)
    ny2, nu2 = sys2.ny, sys2.nu
    nu2 == ny2  || error("The system sys2 must be square")
    ny1, nu1 = sys1.ny, sys1.nu
    nu1 == nu2  || error("The systems sys1 and sys2 must have the same number of inputs")
    nx1 = sys1.nx
    nx2 = sys2.nx
    if nx2 > 0
        A, B, C, D = ssdata([sys1 sys2])
        Ai = [A B[:,1:nu1]; C D[:,1:nu1]]
        Ei = [I zeros(T,nx1+nx2,nu1); zeros(T,nu1,nx1+nx2+nu1)] |> Matrix # TODO: rm call to Matrix when type piracy in https://github.com/JuliaLinearAlgebra/LinearMaps.jl/issues/219 is fixed
        MatrixPencils.isregular(Ai, Ei; atol1, atol2, rtol) || 
            error("The system sys1 is not invertible")
        Bi = [B[:,nu1+1:nu1+nu2]; D[:,nu1+1:nu1+nu2]]
        Ci = [zeros(T,ny1,nx1+nx2) -I] 
        Di = zeros(T,ny1,nu2)
        Ai, Ei, Bi, Ci, Di = MatrixPencils.lsminreal(Ai, Ei, Bi, Ci, Di; fast = true, atol1 = 0, atol2, rtol, contr = true, obs = true, noseig = true)
        if Ei != I
            luE = lu!(Ei, check=false)
            issuccess(luE) || throw(ArgumentError("The system sys1 is not invertible"))
            Ai = luE\Ai
            Bi = luE\Bi
        end
    else
        D1 = T.(sys1.D)
        LUD = lu(D1)
        (norm(D1,Inf) <= atol1 || rcond(LUD.U) <= 10*nu1*eps(real(float(one(T))))) && 
                error("The system sys1 is not invertible")
        Ai, Bi, Ci, Di = ssdata(sys2)
        ldiv!(LUD, Ci); ldiv!(LUD, Di)
    end

    return StateSpace{TE, T}(Ai, Bi, Ci, Di, timeevol) 
end

function /(n::Number, sys::ST) where ST <: AbstractStateSpace
    # Ensure s.D is invertible
    A, B, C, D = ssdata(sys)
    size(D, 1) == size(D, 2) || error("The inverted system must have the same number of inputs and outputs")
    Dinv = try
        inv(D)
    catch
        error("D isn't invertible. If you are trying to form a quotient between two systems `N(s) / D(s)` where the quotient is proper but the inverse of `D(s)` isn't, consider calling `N / D` instead of `N * inv(D)")
    end
    return basetype(ST)(A - B*Dinv*C, B*Dinv, -n*Dinv*C, n*Dinv, sys.timeevol)
end

Base.inv(sys::AbstractStateSpace) = 1/sys
/(sys::ST, n::Number) where ST <: AbstractStateSpace = basetype(ST)(sys.A, sys.B/n, sys.C, sys.D/n, sys.timeevol)
Base.:\(n::Number, sys::ST) where ST <: AbstractStateSpace = basetype(ST)(sys.A, sys.B, sys.C/n, sys.D/n, sys.timeevol)


function Base.adjoint(sys::ST) where ST <: AbstractStateSpace{Continuous}
    return basetype(ST)(-sys.A', -sys.C', sys.B', sys.D') 
end

function Base.adjoint(sys::ST) where ST <: AbstractStateSpace{<:Discrete}
       nx, ny, nu = sys.nx, sys.ny, sys.nu
       T = numeric_type(sys)
       return basetype(ST)(
            [sys.A' sys.C'; zeros(T,ny,nx+ny)], 
            [zeros(T,nx,ny) ; -I], 
            [sys.B' zeros(T,nu,ny)], copy(sys.D'),
            sys.timeevol
        ) 
    end
 end


#####################################################################
##                       Indexing Functions                        ##
#####################################################################
Base.ndims(::AbstractStateSpace) = 2 # NOTE: Also for SISO systems?
Base.size(sys::AbstractStateSpace) = (noutputs(sys), ninputs(sys)) # NOTE: or just size(sys.D)
Base.size(sys::AbstractStateSpace, d::Integer) = d <= 2 ? size(sys)[d] : 1
Base.eltype(::Type{S}) where {S<:AbstractStateSpace} = S
Base.axes(sys::AbstractStateSpace, i::Integer) = Base.OneTo(size(sys, i))

function Base.getindex(sys::ST, inds...) where ST <: AbstractStateSpace
    if size(inds, 1) != 2
        error("Must specify 2 indices to index statespace model")
    end
    rows, cols = index2range(inds...) # FIXME: ControlSystemsBase.index2range(inds...)
    return basetype(ST)(copy(sys.A), sys.B[:, cols], sys.C[rows, :], sys.D[rows, cols], sys.timeevol)
end

function Base.getproperty(sys::AbstractStateSpace, s::Symbol)
    if s === :Ts
        # if !isdiscrete(sys) # NOTE this line seems to be breaking inference of isdiscrete (is there a test for this?)
        if isdiscrete(sys)
            return timeevol(sys).Ts
        else
            @warn "Getting time 0.0 for non-discrete systems is deprecated. Check `isdiscrete` before trying to access time."
            return 0.0
        end
    elseif s === :nx
        return nstates(sys)
    elseif s === :nu
        return ninputs(sys)
    elseif s === :ny
        return noutputs(sys)
    else
        return getfield(sys, s)
    end
end

function Base.propertynames(s::AbstractStateSpace, private::Bool=false)
    (fieldnames(typeof(s))..., :nu, :ny, :nx, (isdiscrete(s) ? (:Ts,) : ())...)
end

#####################################################################
##                        Display Functions                        ##
#####################################################################

_string_mat_with_headers(X) = _string_mat_with_headers(Matrix(X))
function _string_mat_with_headers(X::Matrix)
    #mat = [[""] reshape(cols,1,length(cols));
    #       rows X]
    p = (io, m) -> Base.print_matrix(io, m)
    return replace(sprint(p, X), "\"" => " ")
end

Base.print(io::IO, sys::AbstractStateSpace) = show(io, sys)

function Base.show(io::IO, sys::AbstractStateSpace)
    # Compose the name vectors
    #inputs = format_names(s.inputnames, "u", "?")
    #outputs = format_names(s.outputnames, "y", "?")
    #println(io, "StateSpace:")
    println(io, typeof(sys))
    if nstates(sys) > 0
        #states = format_names(s.statenames, "x", "?")
        println(io, "A = \n", _string_mat_with_headers(sys.A))
        println(io, "B = \n", _string_mat_with_headers(sys.B))
        println(io, "C = \n", _string_mat_with_headers(sys.C))
    end
    println(io, "D = \n", _string_mat_with_headers(sys.D), "\n")
    # Print sample time
    if isdiscrete(sys)
        println(io, "Sample Time: ", sys.Ts, " (seconds)")
    end
    # Print model type
    if iscontinuous(sys)
        print(io, "Continuous-time state-space model")
    else
        print(io, "Discrete-time state-space model")
    end
end




"""
    minreal(sys::T; fast=false, kwargs...)

Minimal realisation algorithm from P. Van Dooreen, The generalized eigenstructure problem in linear system theory, IEEE Transactions on Automatic Control

For information about the options, see `?ControlSystemsBase.MatrixPencils.lsminreal`

See also [`sminreal`](@ref), which is both numerically exact and substantially faster than `minreal`, but with a much more limited potential in removing non-minimal dynamics.
"""
function minreal(sys::T, tol=nothing; fast=false, atol=0.0, kwargs...) where T <: AbstractStateSpace
    A,B,C,D = ssdata(sys)
    if tol !== nothing
        atol == 0 || atol == tol || error("Both positional argument `tol` and keyword argument `atol` were set but were not equal. `tol` is provided for backwards compat and can not be set to another value than `atol`.")
        atol = tol
    end
    Ar, Br, Cr = MatrixPencils.lsminreal(A,B,C; atol, fast, kwargs...)
    if hasfield(T, :sys)
        basetype(T)(ss(Ar,Br,Cr,D), ntuple(i->getfield(sys, i+1), fieldcount(T)-1)...)
    else
        basetype(T)(Ar,Br,Cr,D, ntuple(i->getfield(sys, i+4), fieldcount(T)-4)...)
    end
end


"""
`dsys = diagonalize(s::StateSpace, digits=12)` Diagonalizes the system such that the A-matrix is diagonal.
"""
function diagonalize(s::AbstractStateSpace)
    S,V = eigen(s.A)
    try
        A = diagm(0 => S)
        B = V\s.B
        C = s.C*V
        D = s.D
        return ss(A,B,C,D)
    catch e
        error("System not diagonalizable", e)
    end
end
