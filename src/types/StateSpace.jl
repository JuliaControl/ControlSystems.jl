#####################################################################
##                      Data Type Declarations                     ##
#####################################################################

function state_space_validation(A,B,C,D,Ts::TimeType)
    nx = size(A, 1)
    nu = size(B, 2)
    ny = size(C, 1)

    if size(A, 2) != nx && nx != 0
        error("A has dimentions $(size(A)), but must be square")
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

abstract type AbstractStateSpace <: LTISystem end
struct StateSpace{TimeT<:TimeType, T, MT<:AbstractMatrix{T}} <: AbstractStateSpace
    A::MT
    B::MT
    C::MT
    D::MT
    time::TimeT
    nx::Int
    nu::Int
    ny::Int
    function StateSpace{TimeT, T, MT}(A::MT, B::MT, C::MT, D::MT, Ts::TimeT) where {TimeT<:TimeType, T, MT <: AbstractMatrix{T}}
        nx,nu,ny = state_space_validation(A,B,C,D,Ts)
        new{TimeT, T, MT}(A, B, C, D, Ts, nx, nu, ny)
    end
end
function StateSpace(A::MT, B::MT, C::MT, D::MT, Ts::TimeT) where {TimeT<:TimeType, T, MT <: AbstractMatrix{T}}
    StateSpace{TimeT, T, MT}(A, B, C, D, Ts)
end
# Constructor for Discrete system
function StateSpace(A::MT, B::MT, C::MT, D::MT, Ts::Number) where {T, MT <: AbstractMatrix{T}}
    StateSpace(A, B, C, D, Discrete(Ts))
end
# Constructor for Continuous system
function StateSpace(A::MT, B::MT, C::MT, D::MT) where {T, MT <: AbstractMatrix{T}}
    StateSpace(A, B, C, D, Continuous())
end

""" If D=0 then convert to correct size and type, else, create 1x1 matrix"""
function fix_D_matrix(T::Type,B,C,D)
    if D == 0
        D = fill(zero(T), size(C,1), size(B,2))
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
    C = to_matrix(T, C)
    D = fix_D_matrix(T, B, C, D)
    return A, B, C, D, T
end

const AbstractNumOrArray = Union{Number, AbstractVecOrMat}

# Explicit construction with different types
function StateSpace{TimeT,T,MT}(A::AbstractNumOrArray, B::AbstractNumOrArray, C::AbstractNumOrArray, D::AbstractNumOrArray, Ts::TimeType) where {TimeT, T, MT <: AbstractMatrix{T}}
    D = fix_D_matrix(T, B, C, D)
    return StateSpace{TimeT, T,Matrix{T}}(MT(to_matrix(T, A)), MT(to_matrix(T, B)), MT(to_matrix(T, C)), MT(D), TimeT(Ts))
end

# Explicit conversion
function StateSpace{TimeT,T,MT}(sys::StateSpace) where {TimeT, T, MT <: AbstractMatrix{T}}
    StateSpace{TimeT,T,MT}(sys.A,sys.B,sys.C,sys.D,sys.time)
end

function StateSpace(A::AbstractNumOrArray, B::AbstractNumOrArray, C::AbstractNumOrArray, D::AbstractNumOrArray, Ts::TimeType)
    A, B, C, D, T = to_similar_matrices(A,B,C,D)
    return StateSpace{typeof(Ts),T,Matrix{T}}(A, B, C, D, Ts)
end
# General Discrete constructor
StateSpace(A::AbstractNumOrArray, B::AbstractNumOrArray, C::AbstractNumOrArray, D::AbstractNumOrArray, Ts::Number) =
    StateSpace(A, B, C, D, Discrete(Ts))
# General continuous constructor
StateSpace(A::AbstractNumOrArray, B::AbstractNumOrArray, C::AbstractNumOrArray, D::AbstractNumOrArray) =
    StateSpace(A, B, C, D, Continuous())

# Function for creation of static gain
function StateSpace(D::AbstractArray{T}, Ts::TimeType) where {T<:Number}
    ny, nu = size(D, 1), size(D, 2)
    A = fill(zero(T), 0, 0)
    B = fill(zero(T), 0, nu)
    C = fill(zero(T), ny, 0)
    D = reshape(D, (ny,nu))
    return StateSpace(A, B, C, D, Ts)
end
StateSpace(D::AbstractArray, Ts::Number) = StateSpace(D, Discrete(Ts))
StateSpace(D::AbstractArray) = StateSpace(D, Continuous())

StateSpace(d::Number, Ts::Number; kwargs...) = StateSpace([d], Discrete(Ts))
StateSpace(d::Number; kwargs...) = StateSpace([d], Continuous())


# StateSpace(sys) converts to StateSpace
StateSpace(sys::LTISystem) = convert(StateSpace, sys)

"""
    `sys = ss(A, B, C, D [,Ts])`

Create a state-space model `sys::StateSpace{TimeT, T, MT<:AbstractMatrix{T}}`
where `MT` is the type of matrixes `A,B,C,D`, `T` the element type and TimeT is Continuous or Discrete.

This is a continuous-time model if Ts is omitted.
Otherwise, this is a discrete-time model with sampling period Ts.

`sys = ss(D [, Ts])` specifies a static gain matrix D.
"""
ss(args...;kwargs...) = StateSpace(args...;kwargs...)


struct HeteroStateSpace{TimeT<:TimeType, AT<:AbstractMatrix,BT<:AbstractMatrix,CT<:AbstractMatrix,DT<:AbstractMatrix} <: AbstractStateSpace
    A::AT
    B::BT
    C::CT
    D::DT
    time::TimeT
    nx::Int
    nu::Int
    ny::Int
end
function HeteroStateSpace(A::AT, B::BT,
    C::CT, D::DT, Ts::TimeT) where {TimeT<:TimeType,AT<:AbstractMatrix,BT<:AbstractMatrix,CT<:AbstractMatrix,DT<:AbstractMatrix}
    nx,nu,ny = state_space_validation(A,B,C,D,Ts)
    HeteroStateSpace{TimeT,AT,BT,CT,DT}(A, B, C, D, Ts, nx, nu, ny)
end

# Explicit constructor
function HeteroStateSpace{TimeT,AT,BT,CT,DT}(A, B, C, D, Ts) where {TimeT,AT,BT,CT,DT}
    nx,nu,ny = state_space_validation(A,B,C,D,Ts)
    HeteroStateSpace{TimeT,AT,BT,CT,DT}(AT(A), BT(B), CT(C), DT(D), TimeT(Ts), nx, nu, ny)
end

HeteroStateSpace(s::AbstractStateSpace) = HeteroStateSpace(s.A,s.B,s.C,s.D,s.time)

# Base constructor
function HeteroStateSpace(A::AbstractNumOrArray, B::AbstractNumOrArray, C::AbstractNumOrArray, D::AbstractNumOrArray, Ts::TimeType)
    A = to_abstract_matrix(A)
    B = to_abstract_matrix(B)
    C = to_abstract_matrix(C)
    if D == 0
        D = fill(zero(eltype(C)), size(C,1), size(B,2))
    else
        D = to_abstract_matrix(D)
    end
    return HeteroStateSpace{typeof(Ts),typeof(A),typeof(B),typeof(C),typeof(D)}(A, B, C, D, Ts)
end

HeteroStateSpace(A::AbstractNumOrArray, B::AbstractNumOrArray, C::AbstractNumOrArray, D::AbstractNumOrArray, Ts::Number) =
    HeteroStateSpace(A, B, C, D, Discrete(Ts))
HeteroStateSpace(A::AbstractNumOrArray, B::AbstractNumOrArray, C::AbstractNumOrArray, D::AbstractNumOrArray) =
    HeteroStateSpace(A, B, C, D, Continuous())


# Function for creation of static gain
function HeteroStateSpace(D::AbstractArray{T}, Ts::TimeType) where {T<:Number}
    ny, nu = size(D, 1), size(D, 2)
    A = fill(zero(T), 0, 0)
    B = fill(zero(T), 0, nu)
    C = fill(zero(T), ny, 0)

    return HeteroStateSpace(A, B, C, D, Ts)
end

HeteroStateSpace(D::AbstractArray{T}, Ts::Number) where {T<:Number} = HeteroStateSpace(D, Discrete(Ts))
HeteroStateSpace(D::AbstractArray{T}) where {T<:Number} = HeteroStateSpace(D, Continuous())

HeteroStateSpace(d::Number, Ts::Number; kwargs...) = HeteroStateSpace([d], Discrete(Ts))
HeteroStateSpace(d::Number; kwargs...) = HeteroStateSpace([d], Continuous())

# HeteroStateSpace(sys) converts to HeteroStateSpace
HeteroStateSpace(sys::LTISystem) = convert(HeteroStateSpace, sys)

# Getter functions ###################################################

ssdata(sys::AbstractStateSpace) = sys.A, sys.B, sys.C, sys.D

# Funtions for number of intputs, outputs and states
ninputs(sys::AbstractStateSpace) = size(sys.D, 2)
noutputs(sys::AbstractStateSpace) = size(sys.D, 1)
nstates(sys::AbstractStateSpace) = size(sys.A, 1)

isproper(sys::AbstractStateSpace) = iszero(sys.D)

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
function +(s1::StateSpace{TimeT,T,MT}, s2::StateSpace{TimeT,T,MT}) where {TimeT,T, MT}
    #Ensure systems have same dimensions
    if size(s1) != size(s2)
        error("Systems have different shapes.")
    end
    Ts = common_sample_time(s1.time,s2.time)

    A = [s1.A                   fill(zero(T), nstates(s1), nstates(s2));
         fill(zero(T), nstates(s2), nstates(s1))        s2.A]
    B = [s1.B ; s2.B]
    C = [s1.C s2.C;]
    D = [s1.D + s2.D;]

    return StateSpace{TimeT,T,MT}(A, B, C, D, Ts)
end

function +(s1::HeteroStateSpace, s2::HeteroStateSpace)
    #Ensure systems have same dimensions
    if size(s1) != size(s2)
        error("Systems have different shapes.")
    end
    Ts = common_sample_time(s1.time,s2.time)
    T = promote_type(eltype(s1.A),eltype(s2.A))
    A = [s1.A                   fill(zero(T), nstates(s1), nstates(s2));
         fill(zero(T), nstates(s2), nstates(s1))        s2.A]
    B = [s1.B ; s2.B]
    C = [s1.C s2.C;]
    D = [s1.D + s2.D;]

    return HeteroStateSpace(A, B, C, D, Ts)
end

+(sys::ST, n::Number) where ST <: AbstractStateSpace = ST(sys.A, sys.B, sys.C, sys.D .+ n, sys.time)
+(n::Number, sys::ST) where ST <: AbstractStateSpace = +(sys, n)

## SUBTRACTION ##
-(sys1::AbstractStateSpace, sys2::AbstractStateSpace) = +(sys1, -sys2)
-(sys::AbstractStateSpace, n::Number) = +(sys, -n)
-(n::Number, sys::AbstractStateSpace) = +(-sys, n)

## NEGATION ##
-(sys::ST) where ST <: AbstractStateSpace = ST(sys.A, sys.B, -sys.C, -sys.D, sys.time)

## MULTIPLICATION ##
function *(sys1::StateSpace{TimeT,T,MT}, sys2::StateSpace{TimeT,T,MT}) where {TimeT,T, MT}
    #Check dimension alignment
    #Note: sys1*sys2 = y <- sys1 <- sys2 <- u
    if sys1.nu != sys2.ny
        error("sys1*sys2: sys1 must have same number of inputs as sys2 has outputs")
    end
    Ts = common_sample_time(sys1.time,sys2.time)

    A = [sys1.A    sys1.B*sys2.C;
         fill(zero(T), sys2.nx, sys1.nx)  sys2.A]
    B = [sys1.B*sys2.D ; sys2.B]
    C = [sys1.C   sys1.D*sys2.C;]
    D = [sys1.D*sys2.D;]
    return StateSpace{TimeT,T,MT}(A, B, C, D, Ts)
end

function *(sys1::HeteroStateSpace, sys2::HeteroStateSpace)
    #Check dimension alignment
    #Note: sys1*sys2 = y <- sys1 <- sys2 <- u
    if sys1.nu != sys2.ny
        error("sys1*sys2: sys1 must have same number of inputs as sys2 has outputs")
    end
    Ts = common_sample_time(sys1.time,sys2.time)
    T = promote_type(eltype(sys1.A),eltype(sys2.A))
    A = [sys1.A    sys1.B*sys2.C;
         fill(zero(T), sys2.nx, sys1.nx)  sys2.A]
    B = [sys1.B*sys2.D ; sys2.B]
    C = [sys1.C   sys1.D*sys2.C;]
    D = [sys1.D*sys2.D;]

    return HeteroStateSpace(A, B, C, D, Ts)
end

*(sys::ST, n::Number) where ST <: AbstractStateSpace = StateSpace(sys.A, sys.B, sys.C*n, sys.D*n, sys.time)
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
    return ST(A - B*Dinv*C, B*Dinv, -n*Dinv*C, n*Dinv, sys.time)
end

Base.inv(sys::AbstractStateSpace) = 1/sys
/(sys::ST, n::Number) where ST <: AbstractStateSpace = ST(sys.A, sys.B, sys.C/n, sys.D/n, sys.time)

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
    return ST(copy(sys.A), sys.B[:, cols], sys.C[rows, :], sys.D[rows, cols], sys.time)
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
        println(io, "Sample Time: ", sampletime(sys), " (seconds)")
    end
    # Print model type
    if iscontinuous(sys)
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
