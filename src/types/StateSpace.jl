#####################################################################
##                      Data Type Declarations                     ##
#####################################################################

function state_space_validation(A,B,C,D,Ts)
    nx = size(A, 1)
    nu = size(B, 2)
    ny = size(C, 1)

    if size(A, 2) != nx && nx != 0
        error("A must be square")
    elseif size(B, 1) != nx
        error("B must have the same row size as A")
    elseif size(C, 2) != nx
        error("C must have the same column size as A")
    elseif nu != size(D, 2)
        error("D must have the same column size as B")
    elseif ny != size(D, 1)
        error("D must have the same row size as C")
    end

    # Validate sampling time
    if Ts < 0 && Ts != -1
        error("Ts must be either a positive number, 0
               (continuous system), or -1 (unspecified)")
    end
    nx,nu,ny
end

abstract type AbstractStateSpace <: LTISystem end
struct StateSpace{T, MT<:AbstractMatrix{T}} <: AbstractStateSpace
    A::MT
    B::MT
    C::MT
    D::MT
    Ts::Float64
    nx::Int
    nu::Int
    ny::Int
    function StateSpace{T, MT}(A::MT, B::MT,
            C::MT, D::MT, Ts::Float64) where {T, MT <: AbstractMatrix{T}}
        nx,nu,ny = state_space_validation(A,B,C,D,Ts)
        new{T, MT}(A, B, C, D, Ts, nx, nu, ny)
    end
end

const AbstractNumOrArray = Union{Number, AbstractVecOrMat}

function StateSpace{T,MT}(A::AbstractNumOrArray, B::AbstractNumOrArray, C::AbstractNumOrArray, D::AbstractNumOrArray, Ts::Real) where {T, MT <: AbstractMatrix{T}}
    if D == 0
        D = fill(zero(T), size(C,1), size(B,2))
    else
        D = to_matrix(T, D)
    end
    return StateSpace{T,Matrix{T}}(MT(to_matrix(T, A)), MT(to_matrix(T, B)), MT(to_matrix(T, C)), MT(D), Float64(Ts))
end

function StateSpace{T,MT}(sys::StateSpace) where {T, MT <: AbstractMatrix{T}}
    StateSpace{T,MT}(sys.A,sys.B,sys.C,sys.D,sys.Ts)
end

function StateSpace(A::AbstractNumOrArray, B::AbstractNumOrArray, C::AbstractNumOrArray, D::AbstractNumOrArray, Ts::Real=0)
    T = promote_type(eltype(A),eltype(B),eltype(C),eltype(D))
    A = to_matrix(T, A)
    B = to_matrix(T, B)
    C = to_matrix(T, C)
    if D == 0
        D = fill(zero(T), size(C,1), size(B,2))
    else
        D = to_matrix(T, D)
    end
    return StateSpace{T,Matrix{T}}(A, B, C, D, Float64(Ts))
end


# Function for creation of static gain
function StateSpace(D::AbstractArray{T}, Ts::Real=0) where {T<:Number}
    ny, nu = size(D, 1), size(D, 2)
    A = fill(zero(T), 0, 0)
    B = fill(zero(T), 0, nu)
    C = fill(zero(T), ny, 0)

    return StateSpace(A, B, C, D, Ts)
end
StateSpace(d::Number, Ts::Real=0; kwargs...) = StateSpace([d], Ts)

# StateSpace(sys) converts to StateSpace
StateSpace(sys::LTISystem) = convert(StateSpace, sys)

"""
    `sys = ss(A, B, C, D, Ts=0)`

Create a state-space model `sys::StateSpace{T, MT<:AbstractMatrix{T}}`
where `MT` is the type of matrixes `A,B,C,D` and `T` the element type.

This is a continuous-time model if Ts is omitted or set to 0.
Otherwise, this is a discrete-time model with sampling period Ts.
Set Ts=-1 for a discrete-time model with unspecified sampling period.

`sys = ss(D[, Ts, ...])` specifies a static gain matrix D.
"""
ss(args...;kwargs...) = StateSpace(args...;kwargs...)


struct HeteroStateSpace{AT<:AbstractVecOrMat,BT<:AbstractVecOrMat,CT<:AbstractVecOrMat,DT<:AbstractVecOrMat} <: AbstractStateSpace
    A::AT
    B::BT
    C::CT
    D::DT
    Ts::Float64
    nx::Int
    nu::Int
    ny::Int
end
function HeteroStateSpace(A::AT, B::BT,
    C::CT, D::DT, Ts::Float64=0) where {AT<:AbstractVecOrMat,BT<:AbstractVecOrMat,CT<:AbstractVecOrMat,DT<:AbstractVecOrMat}
    nx,nu,ny = state_space_validation(A,B,C,D,Ts)
    HeteroStateSpace{AT,BT,CT,DT}(A, B, C, D, Ts, nx, nu, ny)
end

function HeteroStateSpace{AT,BT,CT,DT}(A, B, C, D, Ts::Float64=0) where {AT,BT,CT,DT}
    nx,nu,ny = state_space_validation(A,B,C,D,Ts)
    HeteroStateSpace{AT,BT,CT,DT}(AT(A), BT(B), CT(C), DT(D), Ts, nx, nu, ny)
end

HeteroStateSpace(s::AbstractStateSpace) = HeteroStateSpace(s.A,s.B,s.C,s.D,s.Ts)

function HeteroStateSpace(A::AbstractNumOrArray, B::AbstractNumOrArray, C::AbstractNumOrArray, D::AbstractNumOrArray, Ts::Real=0)
    A = to_matrix(eltype(A), A)
    B = to_matrix(eltype(B), B)
    C = to_matrix(eltype(C), C)
    if D == 0
        D = fill(zero(eltype(C)), size(C,1), size(B,2))
    else
        D = to_matrix(eltype(D), D)
    end
    return HeteroStateSpace{typeof(A),typeof(B),typeof(C),typeof(D)}(A, B, C, D, Float64(Ts))
end


# Function for creation of static gain
function HeteroStateSpace(D::AbstractArray{T}, Ts::Real=0) where {T<:Number}
    ny, nu = size(D, 1), size(D, 2)
    A = fill(zero(T), 0, 0)
    B = fill(zero(T), 0, nu)
    C = fill(zero(T), ny, 0)

    return HeteroStateSpace(A, B, C, D, Ts)
end
HeteroStateSpace(d::Number, Ts::Real=0; kwargs...) = HeteroStateSpace([d], Ts)

# HeteroStateSpace(sys) converts to HeteroStateSpace
HeteroStateSpace(sys::LTISystem) = convert(HeteroStateSpace, sys)

# Getter functions ###################################################

ssdata(sys::AbstractStateSpace) = sys.A, sys.B, sys.C, sys.D

# Funtions for number of intputs, outputs and states
ninputs(sys::AbstractStateSpace) = size(sys.D, 2)
noutputs(sys::AbstractStateSpace) = size(sys.D, 1)
nstates(sys::AbstractStateSpace) = size(sys.A, 1)

#####################################################################
##                         Math Operators                          ##
#####################################################################

## EQUALITY ##
function ==(sys1::ST1, sys2::ST2) where {ST1<:AbstractStateSpace,ST2<:AbstractStateSpace}
    fieldnames(ST1) == fieldnames(ST2) || (return false)
    return all(getfield(sys1, f) == getfield(sys2, f) for f in fieldnames(ST1))
end

## Approximate ##
function isapprox(sys1::ST1, sys2::ST2) where {ST1<:AbstractStateSpace,ST2<:AbstractStateSpace}
    fieldnames(ST1) == fieldnames(ST2) || (return false)
    return all(getfield(sys1, f) ≈ getfield(sys2, f) for f in fieldnames(ST1))
end

## ADDITION ##
function +(s1::StateSpace{T,MT}, s2::StateSpace{T,MT}) where {T, MT}
    #Ensure systems have same dimensions
    if size(s1) != size(s2)
        error("Systems have different shapes.")
    elseif s1.Ts != s2.Ts
        error("Sampling time mismatch")
    end

    A = [s1.A                   fill(zero(T), nstates(s1), nstates(s2));
         fill(zero(T), nstates(s2), nstates(s1))        s2.A]
    B = [s1.B ; s2.B]
    C = [s1.C s2.C;]
    D = [s1.D + s2.D;]

    return StateSpace(A, B, C, D, s1.Ts)
end

function +(s1::HeteroStateSpace, s2::HeteroStateSpace)
    #Ensure systems have same dimensions
    if size(s1) != size(s2)
        error("Systems have different shapes.")
    elseif s1.Ts != s2.Ts
        error("Sampling time mismatch")
    end
    T = promote_type(eltype(s1.A),eltype(s2.A))
    A = [s1.A                   fill(zero(T), nstates(s1), nstates(s2));
         fill(zero(T), nstates(s2), nstates(s1))        s2.A]
    B = [s1.B ; s2.B]
    C = [s1.C s2.C;]
    D = [s1.D + s2.D;]

    return HeteroStateSpace(A, B, C, D, s1.Ts)
end

+(sys::ST, n::Number) where ST <: AbstractStateSpace = ST(sys.A, sys.B, sys.C, sys.D .+ n, sys.Ts)
+(n::Number, sys::ST) where ST <: AbstractStateSpace = +(sys, n)

## SUBTRACTION ##
-(sys1::AbstractStateSpace, sys2::AbstractStateSpace) = +(sys1, -sys2)
-(sys::AbstractStateSpace, n::Number) = +(sys, -n)
-(n::Number, sys::AbstractStateSpace) = +(-sys, n)

## NEGATION ##
-(sys::ST) where ST <: AbstractStateSpace = ST(sys.A, sys.B, -sys.C, -sys.D, sys.Ts)

## MULTIPLICATION ##
function *(sys1::StateSpace{T,MT}, sys2::StateSpace{T,MT}) where {T, MT}
    #Check dimension alignment
    #Note: sys1*sys2 = y <- sys1 <- sys2 <- u
    if sys1.nu != sys2.ny
        error("sys1*sys2: sys1 must have same number of inputs as sys2 has outputs")
    elseif sys1.Ts != sys2.Ts
        error("Sampling time mismatch")
    end

    A = [sys1.A    sys1.B*sys2.C;
         fill(zero(T), sys2.nx, sys1.nx)  sys2.A]
    B = [sys1.B*sys2.D ; sys2.B]
    C = [sys1.C   sys1.D*sys2.C;]
    D = [sys1.D*sys2.D;]

    return StateSpace{T,MT}(A, B, C, D, sys2.Ts)
end

function *(sys1::HeteroStateSpace, sys2::HeteroStateSpace)
    #Check dimension alignment
    #Note: sys1*sys2 = y <- sys1 <- sys2 <- u
    if sys1.nu != sys2.ny
        error("sys1*sys2: sys1 must have same number of inputs as sys2 has outputs")
    elseif sys1.Ts != sys2.Ts
        error("Sampling time mismatch")
    end
    T = promote_type(eltype(sys1.A),eltype(sys2.A))
    A = [sys1.A    sys1.B*sys2.C;
         fill(zero(T), sys2.nx, sys1.nx)  sys2.A]
    B = [sys1.B*sys2.D ; sys2.B]
    C = [sys1.C   sys1.D*sys2.C;]
    D = [sys1.D*sys2.D;]

    return HeteroStateSpace(A, B, C, D, sys2.Ts)
end

*(sys::ST, n::Number) where ST <: AbstractStateSpace = StateSpace(sys.A, sys.B, sys.C*n, sys.D*n, sys.Ts)
*(n::Number, sys::AbstractStateSpace) = *(sys, n)

## DIVISION ##
/(sys1::AbstractStateSpace, sys2::AbstractStateSpace) = sys1*inv(sys2)

function /(n::Number, sys::ST) where ST <: AbstractStateSpace
    # Ensure s.D is invertible
    A, B, C, D = ssdata(sys)
    Dinv = try
        inv(D)
    catch
        error("D isn't invertible")
    end
    return ST(A - B*Dinv*C, B*Dinv, -n*Dinv*C, n*Dinv, sys.Ts)
end

Base.inv(sys::AbstractStateSpace) = 1/sys
/(sys::ST, n::Number) where ST <: AbstractStateSpace = ST(sys.A, sys.B, sys.C/n, sys.D/n, sys.Ts)

Base.:^(sys::AbstractStateSpace, p::Integer) = Base.power_by_squaring(sys, p)


#####################################################################
##                       Indexing Functions                        ##
#####################################################################
Base.ndims(::AbstractStateSpace) = 2 # NOTE: Also for SISO systems?
Base.size(sys::AbstractStateSpace) = (noutputs(sys), ninputs(sys)) # NOTE: or just size(sys.D)
Base.size(sys::AbstractStateSpace, d::Integer) = d <= 2 ? size(sys)[d] : 1
Base.eltype(::Type{S}) where {S<:AbstractStateSpace} = S

function Base.getindex(sys::ST, inds...) where ST <: AbstractStateSpace
    if size(inds, 1) != 2
        error("Must specify 2 indices to index statespace model")
    end
    rows, cols = index2range(inds...) # FIXME: ControlSystems.index2range(inds...)
    return ST(copy(sys.A), sys.B[:, cols], sys.C[rows, :], sys.D[rows, cols], sys.Ts)
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
    if sys.Ts > 0
        println(io, "Sample Time: ", sys.Ts, " (seconds)")
    elseif sys.Ts == -1
        println(io, "Sample Time: unspecified")
    end
    # Print model type
    if nstates(sys) == 0
        print(io, "Static gain") # NOTE: Not quite...still has a time type
    elseif iscontinuous(sys)
        print(io, "Continuous-time state-space model")
    else
        print(io, "Discrete-time state-space model")
    end
end




"""
`minsys = minreal(s::StateSpace, tol=sqrt(eps()))` is implemented via `baltrunc` and returns a system on diagonal form.
"""
function minreal(s::AbstractStateSpace, tol=sqrt(eps()))
    s = baltrunc(s, atol=tol, rtol = 0)[1]
    try
        return diagonalize(s)
    catch
        error("Minreal only implemented for diagonalizable systems.")
    end
end



"""
`dsys = diagonalize(s::StateSpace, digits=12)` Diagonalizes the system such that the A-matrix is diagonal.
"""
function diagonalize(s::AbstractStateSpace, digits = 12)
    r = x -> round(x, digits=digits)
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
