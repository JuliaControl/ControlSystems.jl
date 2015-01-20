#####################################################################
##                      Data Type Declarations                     ##
#####################################################################

type StateSpace <: LTISystem
    A::Matrix{Float64}
    B::Matrix{Float64}
    C::Matrix{Float64}
    D::Matrix{Float64}
    Ts::Float64
    nx::Int
    nu::Int
    ny::Int
    statenames::Vector{UTF8String}
    inputnames::Vector{UTF8String}
    outputnames::Vector{UTF8String}

    function StateSpace(A::Matrix{Float64}, B::Matrix{Float64},
            C::Matrix{Float64}, D::Matrix{Float64}, Ts::Float64,
            statenames::Vector{UTF8String}, inputnames::Vector{UTF8String},
            outputnames::Vector{UTF8String})
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

function StateSpace(A::Array, B::Array, C::Array, D::Array, Ts::Real,
        statenames::Vector{UTF8String}, inputnames::Vector{UTF8String},
        outputnames::Vector{UTF8String})
    return StateSpace(float64mat(A), float64mat(B), float64mat(C),
            float64mat(D), float64(Ts), statenames, inputnames, outputnames)
end

#####################################################################
##                      Constructor Functions                      ##
#####################################################################

function ss(A::Array, B::Array, C::Array, D::Array, Ts::Real=0; kwargs...)
    # Check the kwargs for metadata
    nu = size(B, 2)
    ny, nx = size(C, 1, 2)
    kvs = Dict(kwargs)
    statenames = validate_names(kvs, :statenames, nx)
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
ss(sys::LTISystem) = convert(StateSpace, sys)

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
    C = [s1.C s2.C]
    D = [s1.D + s2.D]

    # Naming strategy: If only one sys is named, use that. If the names are the
    # same, use them. If the names conflict, then they are ignored, and the
    # default "" is used.
    statenames = [s1.statenames; s2.statenames]
    if all(s1.inputnames .== "")
        inputnames = s2.inputnames
    elseif all(s2.inputnames .== "") || (s1.inputnames == s2.inputnames)
        inputnames = s1.inputnames
    else
        inputnames = UTF8String["" for i = 1:s1.ny]
    end
    if all(s1.outputnames .== "")
        outputnames = s2.outputnames
    elseif all(s2.outputnames .== "") || (s1.outputnames == s2.outputnames)
        outputnames = s1.outputnames
    else
        outputnames = UTF8String["" for i = 1:s1.nu]
    end
    return StateSpace(A, B, C, D, s1.Ts, statenames, inputnames, outputnames)
end

+(s::StateSpace, n::Real) = StateSpace(s.A, s.B, s.C, s.D .+ n, s.Ts,
        s.statenames, s.inputnames, s.outputnames)
+(n::Real, s::StateSpace) = +(s, n)

## SUBTRACTION ##
-(s1::StateSpace, s2::StateSpace) = +(s1, -s2)
-(s::StateSpace, n::Real) = +(s, -n)
-(n::Real, s::StateSpace) = +(-s, n)

## NEGATION ##
-(s::StateSpace) = StateSpace(s.A, s.B, -s.C, -s.D, s.Ts, s.statenames,
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
    C = [s1.C   s1.D*s2.C]
    D = [s1.D*s2.D]

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
Base.ndims(s::StateSpace) = 2
Base.size(s::StateSpace) = (s.ny, s.nu)
Base.size(s::StateSpace, d) = d <= 2 ? size(s)[d] : 1

function Base.getindex(s::StateSpace, inds...)
    if size(inds, 1) != 2
        error("Must specify 2 indices to index statespace model")
    end
    rows, cols = inds
    return StateSpace([s.A], [s.B[:, cols]], [s.C[rows, :]], [s.D[rows, cols]],
            s.Ts, [s.statenames], [s.inputnames[cols]], [s.outputnames[rows]])
end

#####################################################################
##                        Display Functions                        ##
#####################################################################

# TODO : this is a very hacky way of handling StateSpace printing.
function _string_mat_with_headers(X::Matrix, cols::Vector{UTF8String},
                                rows::Vector{UTF8String})
    mat = [[""] cols';
           rows X]
    p = (io, mat) -> Base.showarray(io, mat, header=false, repr=false)
    return replace(sprint(p, mat), "\"", " ")
end

Base.print(io::IO, s::StateSpace) = show(io, s)

function Base.show(io::IO, s::StateSpace)
    # Compose the name vectors
    inputs = format_names(s.inputnames, "u", "?")
    outputs = format_names(s.outputnames, "y", "?")
    println(io, "StateSpace:")
    if s.nx > 0
        states = format_names(s.statenames, "x", "?")
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
