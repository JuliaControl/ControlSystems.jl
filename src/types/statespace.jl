#####################################################################
##                      Data Type Declarations                     ##
#####################################################################

type StateSpace{MT} <: LTISystem where MT <: AbstractNumberMatrix
    A::MT
    B::MT
    C::MT
    D::MT
    Ts::Float64
    nx::Int
    nu::Int
    ny::Int

    function StateSpace{MT}(A::MT, B::MT, C::MT, D::MT, Ts::Float64) where MT <: AbstractNumberMatrix
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
        new{MT}(A, B, C, D, Ts, nx, nu, ny)
    end
end

function StateSpace(A::AbstractArray, B::AbstractArray, C::AbstractArray, D::AbstractArray, Ts::Real)
        T = promote_type(eltype(A),eltype(B),eltype(C),eltype(D))
        @assert (typeof(to_matrix(T, A)) == typeof(to_matrix(T, B)) == typeof(to_matrix(T, C)) == typeof(to_matrix(T, D)))
        return StateSpace{Matrix{T}}(to_matrix(T, A), to_matrix(T, B), to_matrix(T, C), to_matrix(T, D), Float64(Ts))
end

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
function ss(A::Array, B::Array, C::Array, D::Array, Ts::Real=0; kwargs...)
    # Check the kwargs for metadata
    nu = size(B, 2)
    ny, nx = size(C, 1, 2)
    # kvs = Dict(kwargs)
    # statenames = validate_names(kvs, :statenames, nx)
    # inputnames = validate_names(kvs, :inputnames, nu)
    # outputnames = validate_names(kvs, :outputnames, ny)
    return StateSpace(A, B, C, D, Ts)
end

# Function for accepting scalars
function ss(A::Union{Real,Array}, B::Union{Real,Array}, C::Union{Real,Array}, D::Union{Real,Array}, args...; kwargs...)
    T = promote_type(eltype(A),eltype(B),eltype(C),eltype(D))
    A = to_matrix(T, A)
    B = to_matrix(T, B)
    C = to_matrix(T, C)
    if D == 0
        D = zeros(T, size(C,1),size(B,2))
    else
        D = to_matrix(T, D)
    end
    ss(A, B, C, D, args..., kwargs...)
end

# Function for creation of static gain
function ss(D::Array, Ts::Real=0; kwargs...)
    ny, nu = size(D, 1, 2)
    A = zeros(eltype(D), 0, 0)
    B = zeros(eltype(D), 0, nu)
    C = zeros(eltype(D), ny, 0)

    return ss(A, B, C, D, Ts, kwargs...)
end
ss(d::Real, Ts::Real=0; kwargs...) = ss([d], Ts, kwargs...)

# ss(sys) converts to StateSpace
ss(sys::LTISystem) = convert(StateSpace, sys)

# Create a random statespace system
function rss(nx::Int, nu::Int=1, ny::Int=1, feedthrough::Bool=true)
    Q = randn(nx, nx)
    A = Q*diagm(-100*abs.(randn(nx)))*Q'
    B = randn(nx, nu)
    C = randn(ny, nx)
    if feedthrough
        D = randn(ny, nu)
    else
        D = zeros(ny, nu)
    end
    return ss(A, B, C, D)
end

#####################################################################
##                         Math Operators                          ##
#####################################################################

## EQUALITY ##
function ==(s1::StateSpace, s2::StateSpace)
    fields = [:Ts, :nx, :ny, :nu, :A, :B, :C, :D]
    for field in fields
        if getfield(s1, field) != getfield(s2, field)
            return false
        end
    end
    return true
end

## Approximate ##
function isapprox(s1::StateSpace, s2::StateSpace)
    fieldsApprox = [:Ts, :nx, :ny, :nu, :A, :B, :C, :D]
    for field in fieldsApprox
        if !(getfield(s1, field) ≈ getfield(s2, field))
            return false
        end
    end
    return true
end

## ADDITION ##
function +(s1::StateSpace, s2::StateSpace)
    #Ensure systems have same dimensions
    if size(s1) != size(s2)
        error("Systems have different shapes.")
    elseif s1.Ts != s2.Ts
        error("Sampling time mismatch")
    end

    A = [s1.A zeros(s1.nx, s2.nx);
         zeros(s2.nx, s1.nx) s2.A]
    B = [s1.B ; s2.B]
    C = [s1.C s2.C;]
    D = [s1.D + s2.D;]

    return StateSpace(A, B, C, D, s1.Ts)
end

+(s::StateSpace, n::Real) = StateSpace(s.A, s.B, s.C, s.D .+ n, s.Ts)
+(n::Real, s::StateSpace) = +(s, n)

## SUBTRACTION ##
-(s1::StateSpace, s2::StateSpace) = +(s1, -s2)
-(s::StateSpace, n::Real) = +(s, -n)
-(n::Real, s::StateSpace) = +(-s, n)

## NEGATION ##
-(s::StateSpace) = StateSpace(s.A, s.B, -s.C, -s.D, s.Ts)

## MULTIPLICATION ##
function *(s1::StateSpace, s2::StateSpace)
    #Check dimension alignment
    #Note: s1*s2 = y <- s1 <- s2 <- u
    if s1.nu != s2.ny
        error("s1*s2: s1 must have same number of inputs as s2 has outputs")
    elseif s1.Ts != s2.Ts
        error("Sampling time mismatch")
    end

    A = [s1.A    s1.B*s2.C;
         zeros(s2.nx, s1.nx)  s2.A]
    B = [s1.B*s2.D ; s2.B]
    C = [s1.C   s1.D*s2.C;]
    D = [s1.D*s2.D;]

    return StateSpace(A, B, C, D, s1.Ts)
end

*(s::StateSpace, n::Real) = StateSpace(s.A, s.B, s.C*n, s.D*n, s.Ts)
*(n::Real, s::StateSpace) = *(s, n)

## DIVISION ##
/(s1::StateSpace, s2::StateSpace) = s1*inv(s2)

function /(n::Real, s::StateSpace)
    # Ensure s.D is invertible
    Dinv = try
        inv(s.D)
    catch
        error("D isn't invertible")
    end
    return StateSpace(s.A - s.B*Dinv*s.C, s.B*Dinv, -n*Dinv*s.C, n*Dinv, s.Ts)
end

Base.inv(s::StateSpace) = 1/s
/(s::StateSpace, n::Real) = StateSpace(s.A, s.B, s.C/n, s.D/n, s.Ts)

#####################################################################
##                       Indexing Functions                        ##
#####################################################################
Base.ndims(::StateSpace) = 2
Base.size(s::StateSpace) = (s.ny, s.nu)
Base.size(s::StateSpace, d) = d <= 2 ? size(s)[d] : 1

function Base.getindex(s::StateSpace, inds...)
    if size(inds, 1) != 2
        error("Must specify 2 indices to index statespace model")
    end
    rows, cols = ControlSystems.index2range(inds...)
    return StateSpace([s.A;], [s.B[:, cols];], [s.C[rows, :];], [s.D[rows, cols];], s.Ts)
end

#####################################################################
##                        Display Functions                        ##
#####################################################################

# TODO : this is a very hacky way of handling StateSpace printing.
function _string_mat_with_headers(X::Matrix, cols::Vector{String},
                                rows::Vector{String})
    mat = [[""] reshape(cols,1,length(cols));
           rows X]
    p = (io, m) -> Base.showarray(io, m, false, header=false)
    return replace(sprint(p, mat), "\"", " ")
end

Base.print(io::IO, s::StateSpace) = show(io, s)

function Base.show(io::IO, s::Union{T,NamedSystem{T}}) where {T<:StateSpace}
    # Compose the name vectors
    inputs = format_names(inputnames(s), "u", "?")
    outputs = format_names(outputnames(s), "y", "?")
    println(io, "StateSpace:")
    if s.nx > 0
        states = format_names(fill("", s.nx), "x", "?")
        println(io, "A = \n", _string_mat_with_headers(s.A, states, states))
        println(io, "B = \n", _string_mat_with_headers(s.B, inputs, states))
        println(io, "C = \n", _string_mat_with_headers(s.C, states, outputs))
    end
    println(io, "D = \n", _string_mat_with_headers(s.D, inputs, outputs), "\n")
    # Print sample time
    if s.Ts > 0
        println(io, "Sample Time: ", s.Ts, " (seconds)")
    elseif s.Ts == -1
        println(io, "Sample Time: unspecified")
    end
    # Print model type
    if s.nx == 0
        print(io, "Static gain")
    elseif iscontinuous(s)
        print(io, "Continuous-time state-space model")
    else
        print(io, "Discrete-time state-space model")
    end
end



#####################################################################
##                        Other  Functions                         ##
#####################################################################

"""
`minsys = minreal(s::StateSpace, tol=sqrt(eps()))` is implemented via `baltrunc` and returns a system on diagonal form.
"""
function minreal(s::StateSpace, tol=sqrt(eps()))
    s = baltrunc(s, atol=tol, rtol = 0)[1]
    try
        return diagonalize(s)
    catch
        error("Minreal only implemented for diagonalizable systems.")
    end
end

"""
`dsys = diagonalize(s::StateSpace, digits=12)` Diagonalizes the system such that the A-matrix is diagonal. The result is rounded to `digits` decimal points.
"""
function diagonalize(s::StateSpace, digits = 12)
    r = x -> round(x,digits)
    S,V = eig(s.A)
    try
        A = V\s.A*V     .|> r
        B = V\s.B       .|> r
        C = s.C*V       .|> r
        D = s.D         .|> r
        return ss(A,B,C,D)
    catch e
        error("System not diagonalizable", e)
    end
end
