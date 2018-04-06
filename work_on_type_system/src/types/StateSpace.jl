#####################################################################
##                      Data Type Declarations                     ##
#####################################################################

type StateSpace{T, MT<:AbstractMatrix{T}} <: LTISystem
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
        new{T, MT}(A, B, C, D, Ts, nx, nu, ny)
    end
end

function StateSpace(A::AbstractArray, B::AbstractArray, C::AbstractArray, D::AbstractArray, Ts::Real)
        T = promote_type(eltype(A),eltype(B),eltype(C),eltype(D))
        @assert (typeof(to_matrix(T, A)) == typeof(to_matrix(T, B)) == typeof(to_matrix(T, C)) == typeof(to_matrix(T, D)))
        return StateSpace{T,Matrix{T}}(to_matrix(T, A), to_matrix(T, B), to_matrix(T, C),
            to_matrix(T, D), Float64(Ts))
end

# Getter functions
get_A(sys::StateSpace) = sys.A
get_B(sys::StateSpace) = sys.B
get_C(sys::StateSpace) = sys.C
get_D(sys::StateSpace) = sys.D

get_Ts(sys::StateSpace) = sys.Ts

get_ssdata(sys::StateSpace) = get_A(sys), get_B(sys), get_C(sys), get_D(sys)

# Funtions for number of intputs, outputs and states
ninputs(sys::StateSpace) = size(get_D(sys), 2)
noutputs(sys::StateSpace) = size(get_D(sys), 1)
nstates(sys::StateSpace) = size(get_A(sys), 1)

#####################################################################
##                      Constructor Functions                      ##
#####################################################################

@doc """`ss(A,B,C,D[, Ts, statenames=..., inputnames=..., outputnames=...]) -> sys`

Create a state-space model.
This is a continuous-time model if Ts is omitted or set to 0.
Otherwise, this is a discrete-time model with sampling period Ts.
Set Ts=-1 for a discrete-time model with unspecified sampling period.

State, input and output names: each can be either a vector of strings (one string per dimension),
or a single string (e.g., "x"). In the latter case, an index is automatically appended to identify
the coordinates for each dimension (e.g. "x1", "x2", ...).

`sys = ss(D[, Ts, ...])` specifies a static gain matrix D.""" ->
function ss(A::Array, B::Array, C::Array, D::Array, Ts::Real=0)
    # Check the kwargs for metadata
    nu = size(B, 2)
    ny, nx = size(C, 1, 2)
    return StateSpace(A, B, C, D, Ts)
end

# Function for accepting scalars
function ss(A::Union{Real,Array}, B::Union{Real,Array}, C::Union{Real,Array}, D::Union{Real,Array}, Ts::Real=0)
    T = promote_type(eltype(A),eltype(B),eltype(C),eltype(D))
    A = to_matrix(T, A)
    B = to_matrix(T, B)
    C = to_matrix(T, C)
    if D == 0
        D = fill(zero(T), size(C,1), size(B,2))
    else
        D = to_matrix(T, D)
    end
    ss(A, B, C, D, Ts)
end

# Function for creation of static gain
function ss(D::Array{T}, Ts::Real=0) where {T<:Number}
    ny, nu = size(D, 1, 2)
    A = fill(zero(T), 0, 0)
    B = fill(zero(T), 0, nu)
    C = fill(zero(T), ny, 0)

    return ss(A, B, C, D, Ts)
end
ss(d::Real, Ts::Real=0; kwargs...) = ss([d], Ts)

# ss(sys) converts to StateSpace
ss(sys::LTISystem) = convert(StateSpace, sys)

#####################################################################
##                         Math Operators                          ##
#####################################################################

## EQUALITY ##
function ==(s1::StateSpace, s2::StateSpace)
    return all(getfield(s1, f) == getfield(s2, f) for f in fieldnames(StateSpace))
end

## Approximate ##
function isapprox(s1::StateSpace, s2::StateSpace)
    fieldsApprox = [:Ts, :nx, :ny, :nu, :A, :B, :C, :D]
    fieldsEqual = [:inputnames, :outputnames, :statenames]
    return all(s1.f ≈ s2.f for f in fieldnames(fieldsApprox)) &&
           all(s1.f == s2.f for f in fieldnames(fieldsApprox))
end

## ADDITION ##
function +(s1::StateSpace{T,MT}, s2::StateSpace{T,MT}) where {T, MT}
    #Ensure systems have same dimensions
    if size(s1) != size(s2)
        error("Systems have different shapes.")
    elseif s1.Ts != s2.Ts
        error("Sampling time mismatch")
    end

    A = [s1.A                   fill(zero(T), s1.nx, s2.nx);
         fill(zero(T), s2.nx, s1.nx)        s2.A]
    B = [s1.B ; s2.B]
    C = [s1.C s2.C;]
    D = [s1.D + s2.D;]

    return StateSpace(A, B, C, D, s1.Ts)
end

+(sys::StateSpace, n::Real) = StateSpace(sys.A, sys.B, sys.C, sys.D .+ n, sys.Ts)
+(n::Real, sys::StateSpace) = +(s, n)

## SUBTRACTION ##
-(sys1::StateSpace, sys2::StateSpace) = +(sys1, -sys2)
-(sys::StateSpace, n::Real) = +(sys, -n)
-(n::Real, sys::StateSpace) = +(-sys, n)

## NEGATION ##
-(sys::StateSpace) = StateSpace(sys.A, sys.B, -sys.C, -sys.D, sys.Ts)

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

*(sys::StateSpace, n::Real) = StateSpace(sys.A, sys.B, sys.C*n, sys.D*n, sys.Ts)
*(n::Real, sys::StateSpace) = *(sys, n)

## DIVISION ##
/(sys1::StateSpace, sys2::StateSpace) = sys1*inv(sys2)

function /(n::Real, sys::StateSpace)
    # Ensure s.D is invertible
    A, B, C, D = get_ssdata(sys)
    Dinv = try
        inv(D)
    catch
        error("D isn't invertible")
    end
    return StateSpace(A - B*Dinv*C, B*Dinv, -n*Dinv*C, n*Dinv, get_Ts(sys))
end

Base.inv(sys::StateSpace) = 1/sys
/(sys::StateSpace, n::Real) = StateSpace(sys.A, sys.B, sys.C/n, sys.D/n, sys.Ts)

#####################################################################
##                       Indexing Functions                        ##
#####################################################################
Base.ndims(::StateSpace) = 2 # NOTE: Also for SISO systems?
Base.size(sys::StateSpace) = (noutputs(sys), ninputs(sys)) # NOTE: or just size(get_D(sys))
Base.size(sys::StateSpace, d) = d <= 2 ? size(sys)[d] : 1

function Base.getindex(sys::StateSpace, inds...)
    if size(inds, 1) != 2
        error("Must specify 2 indices to index statespace model")
    end
    rows, cols = index2range(inds...) # FIXME: ControlSystems.index2range(inds...)
    return StateSpace([sys.A;], [sys.B[:, cols];], [sys.C[rows, :];], [sys.D[rows, cols];], sys.Ts)
end

#####################################################################
##                        Display Functions                        ##
#####################################################################

# TODO : this is a very hacky way of handling StateSpace printing.
# function _string_mat_with_headers(X::Matrix, cols::Vector{String},
#                                 rows::Vector{String})
#     mat = [[""] reshape(cols,1,length(cols));
#            rows X]
#     p = (io, m) -> Base.showarray(io, m, false, header=false)
#     return replace(sprint(p, mat), "\"", " ")
# end

function _string_mat_with_headers(X::Matrix)
    #mat = [[""] reshape(cols,1,length(cols));
    #       rows X]
    p = (io, m) -> Base.showarray(io, m, false, header=false)
    return replace(sprint(p, X), "\"", " ")
end

Base.print(io::IO, sys::StateSpace) = show(io, sys)

function Base.show(io::IO, sys::StateSpace)
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
