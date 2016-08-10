#####################################################################
##                      Data Type Declarations                     ##
#####################################################################

type StateSpace{T} <: LTISystem
    A::Matrix{T}
    B::Matrix{T}
    C::Matrix{T}
    D::Matrix{T}
    Ts::Float64
    nx::Int
    nu::Int
    ny::Int
    statenames::Vector{String}
    inputnames::Vector{String}
    outputnames::Vector{String}

    function StateSpace(A::Matrix{T}, B::Matrix{T},
            C::Matrix{T}, D::Matrix{T}, Ts::Float64,
            statenames::Vector{String}, inputnames::Vector{String},
            outputnames::Vector{String})
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
        # Validate names of state, input, and output
        if size(statenames, 1) != nx
            error("Must have same number of statenames as states")
        elseif size(inputnames, 1) != nu
            error("Must have same number of inputnames as inputs")
        elseif size(outputnames, 1) != ny
            error("Must have same number of outputnames as outputs")
        end
        # Validate sampling time
        if Ts < 0 && Ts != -1
            error("Ts must be either a positive number, 0
                   (continuous system), or -1 (unspecified)")
        end
        new(A, B, C, D, Ts, nx, nu, ny, statenames, inputnames, outputnames)
    end
end

# Or is it already known when the constructor is called what type fo
function StateSpace(A::Union{Number,Array}, B::Union{Number,Array}, C::Union{Number,Array}, D::Union{Number,Array}, Ts::Real,
        statenames::Vector{String}, inputnames::Vector{String},
        outputnames::Vector{String})

        (A, B, C, D) = promote(ensure_matrix(A), ensure_matrix(B),
                               ensure_matrix(C), ensure_matrix(D))

        # Why do we need to specify T, isn't the inner constructor more specific anyway?
        T = promote_type(typeof(A[1]), typeof(B[1]), typeof(C[1]), typeof(D[1]))

        return StateSpace{T}(A, B, C, D, map(Float64, Ts),
                        statenames, inputnames, outputnames)
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
function ss(A::Union{Number,Array}, B::Union{Number,Array}, C::Union{Number,Array}, D::Union{Number,Array}, Ts::Real=0; kwargs...)
    # For accepting scalar D
    if D == 0
        D = zeros(size(C,1),size(B,2))
    end

    nu = size(B, 2)
    ny, nx = size(C, 1, 2)
    kvs = Dict(kwargs)
    statenames = validate_names(kvs, :statenames, nx) # Are these needed?, a transfer function shouldn't have to carry along variable names
    inputnames = validate_names(kvs, :inputnames, nu)
    outputnames = validate_names(kvs, :outputnames, ny)

    return StateSpace(A, B, C, D, Ts, statenames, inputnames, outputnames)
end

# Function for creation of static gain
function ss(D::Array, Ts::Real=0; kwargs...)
    ny, nu = size(D, 1, 2)
    A = zeros(0, 0)
    B = zeros(0, nu)
    C = zeros(ny, 0)
    return ss(A, B, C, D, Ts, kwargs...)
end
ss(d::Real, Ts::Real=0; kwargs...) = ss([d], Ts, kwargs...)

# ss(sys) converts to StateSpace
ss(sys::LTISystem) = convert(StateSpace{Float64}, sys) # TODO: Just defaulting to Float64

# Create a random statespace system
function rss(nx::Int, nu::Int=1, ny::Int=1, feedthrough::Bool=true)
    Q = randn(nx, nx)
    A = Q*diagm(-100*abs(randn(nx)))*Q'
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
##                         Conversion                              ##
#####################################################################

# Convert state-space system of one type to another
function Base.convert{T}(::Type{StateSpace{T}}, sys::StateSpace)
    A = convert(Array{T,2}, sys.A)
    B = convert(Array{T,2}, sys.B)
    C = convert(Array{T,2}, sys.C)
    D = convert(Array{T,2}, sys.D)
    return StateSpace(A, B, C, D, sys.Ts, sys.statenames, sys.inputnames, sys.outputnames)
end

#####################################################################
##                         Math Operators                          ##
#####################################################################

## EQUALITY ##
function ==(s1::StateSpace, s2::StateSpace)
    fields = [:Ts, :nx, :ny, :nu, :A, :B, :C, :D, :inputnames, :outputnames,
            :statenames]
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
    fieldsEqual = [:inputnames, :outputnames, :statenames]
    for field in fieldsApprox
        if !(getfield(s1, field) ≈ getfield(s2, field))
            return false
        end
    end
    for field in fieldsEqual
        if getfield(s1, field) != getfield(s2, field)
            return false
        end
    end
    return true
end

## ADDITION ##
function +{T}(s1::StateSpace{T}, s2::StateSpace{T})
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

    # Naming strategy: If only one sys is named, use that. If the names are the
    # same, use them. If the names conflict, then they are ignored, and the
    # default "" is used.
    statenames = [s1.statenames; s2.statenames]
    if all(s1.inputnames .== "")
        inputnames = s2.inputnames
    elseif all(s2.inputnames .== "") || (s1.inputnames == s2.inputnames)
        inputnames = s1.inputnames
    else
        inputnames = fill(String(""),s1.ny)
    end
    if all(s1.outputnames .== "")
        outputnames = s2.outputnames
    elseif all(s2.outputnames .== "") || (s1.outputnames == s2.outputnames)
        outputnames = s1.outputnames
    else
        outputnames = fill(String(""),s1.nu)
    end
    return StateSpace(A, B, C, D, s1.Ts, statenames, inputnames, outputnames)
end

+(s::StateSpace, n::Number) = StateSpace(s.A, s.B, s.C, s.D .+ n, s.Ts,
        s.statenames, s.inputnames, s.outputnames)
+(n::Real, s::StateSpace) = +(s, n)

## SUBTRACTION ##
-(s1::StateSpace, s2::StateSpace) = +(s1, -s2)
-(s::StateSpace, n::Number) = +(s, -n)
-(n::Number, s::StateSpace) = +(-s, n)

## NEGATION ##
-{T}(s::StateSpace{T}) = StateSpace{T}(s.A, s.B, -s.C, -s.D, s.Ts, s.statenames,
        s.inputnames, s.outputnames)

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

    statenames = [s1.statenames; s2.statenames]
    return StateSpace(A, B, C, D, s1.Ts, statenames, s2.inputnames,
            s1.outputnames)
end

*(s::StateSpace, n::Real) = StateSpace(s.A, s.B, s.C*n, s.D*n, s.Ts,
        s.statenames, s.inputnames, s.outputnames)
*(n::Real, s::StateSpace) = *(s, n)

## DIVISION ##
/(s1::StateSpace, s2::StateSpace) = s1*(1/s2)

function /(n::Real, s::StateSpace)
    # Ensure s.D is invertible
    Dinv = try
        inv(s.D)
    catch
        error("D isn't invertible")
    end
    return StateSpace(s.A - s.B*Dinv*s.C, s.B*Dinv, -n*Dinv*s.C, n*Dinv, s.Ts,
            s.statenames, s.outputnames, s.inputnames)
end

/(s::StateSpace, n::Real) = StateSpace(s.A, s.B, s.C/n, s.D/n, s.Ts,
        s.statenames, s.inputnames, s.outputnames)

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
    return StateSpace([s.A;], [s.B[:, cols];], [s.C[rows, :];], [s.D[rows, cols];],
            s.Ts, [s.statenames;], [s.inputnames[cols];], [s.outputnames[rows];])
end

#####################################################################
##                        Display Functions                        ##
#####################################################################

# TODO : this is a very hacky way of handling StateSpace printing.
function _string_mat_with_headers(X::Matrix, col_labels::Vector{String},
                                row_labels::Vector{String})
    mat = ["" reshape(col_labels,1,length(col_labels));
           row_labels X]
    p = (io, m) -> Base.showarray(io, m, false, header=false)
    return replace(sprint(p, mat), "\"", " ")
end

Base.print(io::IO, s::StateSpace) = show(io, s)

function Base.show(io::IO, sys::StateSpace)
    # Compose the name vectors
    inputnames = format_names(sys.inputnames, "u", "?")
    outputnames = format_names(sys.outputnames, "y", "?")
    println(io, "StateSpace{", typeof(sys.A[1]), "}:")
    if sys.nx > 0
        statenames = format_names(sys.statenames, "x", "?")
        println(io, "A = \n", _string_mat_with_headers(sys.A, statenames, statenames))
        println(io, "B = \n", _string_mat_with_headers(sys.B, inputnames, statenames))
        println(io, "C = \n", _string_mat_with_headers(sys.C, statenames, outputnames))
    end
    println(io, "D = \n", _string_mat_with_headers(sys.D, inputnames, outputnames), "\n")
    # Print sample time
    if sys.Ts > 0
        println(io, "Sample Time: ", sys.Ts, " (seconds)")
    elseif sys.Ts == -1
        println(io, "Sample Time: unspecified")
    end
    # Print model type
    if sys.nx == 0
        print(io, "Static gain")
    elseif iscontinuous(sys)
        print(io, "Continuous-time state-space model")
    else
        print(io, "Discrete-time state-space model")
    end
end
