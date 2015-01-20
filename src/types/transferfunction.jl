include("polys.jl")
include("sisotf.jl")
#####################################################################
##                      Data Type Declarations                     ##
#####################################################################

type TransferFunction <: LTISystem
    matrix::Matrix{SisoTf}
    Ts::Float64
    nu::Int
    ny::Int
    inputnames::Vector{UTF8String}
    outputnames::Vector{UTF8String}
    function TransferFunction(matrix::Matrix{SisoTf}, Ts::Float64,
            inputnames::Vector{UTF8String}, outputnames::Vector{UTF8String})
        # Validate size of input and output names
        ny, nu = size(matrix)
        if size(inputnames, 1) != nu
            error("Must have same number of inputnames as inputs")
        elseif size(outputnames, 1) != ny
            error("Must have same number of outputnames as outputs")
        end
        # Validate sampling time
        if Ts < 0 && Ts != -1
            error("Ts must be either a positive number, 0
                (continuous system), or -1 (unspecified)")
        end
        return new(matrix, Ts, nu, ny, inputnames, outputnames)
    end
end

#####################################################################
##                      Constructor Functions                      ##
#####################################################################

function tf(num::VecOrMat{Vector}, den::VecOrMat{Vector}, Ts::Real=0; kwargs...)
    # Validate input and output dimensions match
    ny, nu = size(num, 1, 2)
    if (ny, nu) != size(den, 1, 2)
        error("num and den dimensions must match")
    end
    matrix = Array(SisoTf, ny, nu)
    for o=1:ny
        for i=1:nu
            matrix[o, i] = SisoTf(num[o, i], den[o, i])
        end
    end
    kvs = Dict(kwargs)
    inputnames = validate_names(kvs, :inputnames, nu)
    outputnames = validate_names(kvs, :outputnames, ny)
    return TransferFunction(matrix, float64(Ts), inputnames, outputnames)
end

tf(num::Vector, den::Vector, args...) =
tf(reshape(Vector[num], 1, 1), reshape(Vector[den], 1, 1), args...)

#####################################################################
##                          Misc. Functions                        ##
#####################################################################

## INDEXING ##
Base.ndims(::TransferFunction) = 2
Base.size(t::TransferFunction) = (t.ny, t.nu)
Base.size(t::TransferFunction, d) = d <= 2 ? size(t)[d] : 1

function Base.getindex(t::TransferFunction, inds...)
    if size(inds, 1) != 2
        error("Must specify 2 indices to index TransferFunction model")
    end
    rows, cols = inds
    mat = Array(SisoTf, length(rows), length(cols))
    mat[:, :] = t.matrix[rows, cols]
    return TransferFunction(mat, t.Ts, [t.inputnames] [t.outputnames])
end

function Base.copy(t::TransferFunction)
    matrix = copy(t.matrix)
    inputnames = copy(t.inputnames)
    outputnames = copy(t.outputnames)
    return TransferFunction(matrix, t.Ts, inputnames, outputnames)
end

#####################################################################
##                         Math Operators                          ##
#####################################################################

## EQUALITY ##
function ==(t1::TransferFunction, t2::TransferFunction)
    fields = [:Ts, :ny, :nu, :inputnames, :outputnames, :statenames, :matrix]
    for field in fields
        if getfield(t1, field) != getfield(t2, field)
            return false
        end
    end
    return true
end

## ADDITION ##
function +(t1::TransferFunction, t2::TransferFunction)
    if size(t1) != size(t2)
        error("Systems have different shapes.")
    elseif t1.Ts != t2.Ts
        error("Sampling time mismatch")
    end

    # Naming strategy: If only one sys is named, use that. If the names are the
    # same, use them. If the names conflict, then they are ignored, and the
    # default "" is used.
    if all(t1.inputnames .== "")
        inputnames = t2.inputnames
    elseif all(t2.inputnames .== "") || (t1.inputnames == t2.inputnames)
        inputnames = t1.inputnames
    else
        inputnames = UTF8String["" for i = 1:t1.ny]
    end
    if all(t1.outputnames .== "")
        outputnames = t2.outputnames
    elseif all(t2.outputnames .== "") || (t1.outputnames == t2.outputnames)
        outputnames = t1.outputnames
    else
        outputnames = UTF8String["" for i = 1:t1.nu]
    end
    t1, t2 = promote(t1, t2)
    matrix = t1.matrix + t2.matrix
    return TransferFunction(matrix, t1.Ts, inputnames, outputnames)
end

+(t::TransferFunction, n::Real) = TransferFunction(t.matrix + n, t.Ts,
        t.inputnames, t.outputnames)
+(n::Real, t::TransferFunction) = +(t, n)

## SUBTRACTION ##
-(n::Real, t::TransferFunction) = TransferFunction(n - t.matrix, t.Ts,
        t.inputnames, t.outputnames)
-(t1::TransferFunction, t2::TransferFunction) = +(t1, -t2)
-(t::TransferFunction, n::Real) = +(t, -n)

## NEGATION ##
-(t::TransferFunction) = TransferFunction(-t.matrix, t.Ts, t.inputnames,
        t.outputnames)

## MULTIPLICATION ##
function *(t1::TransferFunction, t2::TransferFunction)
    # Note: t1*t2 = y <- t1 <- t2 <- u
    if t1.nu != t2.ny
        error("t1*t2: t1 must have same number of inputs as t2 has outputs")
    elseif t1.Ts != t2.Ts
        error("Sampling time mismatch")
    end
    matrix = t1.matrix*t2.matrix
    return TransferFunction(matrix, t1.Ts, t2.inputnames, t1.outputnames)
end

*(t::TransferFunction, n::Real) = TransferFunction(n*t.matrix, t.Ts,
        t.inputnames, t.outputnames)
*(n::Real, t::TransferFunction) = *(t, n)

## DIVISION ##
function /(n::Real, t::TransferFunction)
    if issiso(t)
        matrix = reshape([n/t.matrix[1,1]], 1, 1)
    else
        error("MIMO TransferFunction inversion isn't implemented yet")
    end
    return TransferFunction(matrix, t.Ts, t.outputnames, t.inputnames)
end
/(t::TransferFunction, n::Real) = t*(1/n)
/(t1::TransferFunction, t2::TransferFunction) = t1*(1/t2)

#####################################################################
##                        Display Functions                        ##
#####################################################################

Base.print(io::IO, t::TransferFunction) = show(io, t)

function Base.show(io::IO, t::TransferFunction)
    # Compose the name vectors
    inputs = format_names(t.inputnames, "Input ", "?")
    outputs = format_names(t.outputnames, "Output ", "?")
    println(io, "TransferFunction:")
    var = iscontinuous(t) ? :s : :z
    for i=1:t.nu
        for o=1:t.ny
            if !issiso(t)
                println(io, inputs[i], " to ", outputs[o])
            end
            print_sisotf(io, t.matrix[o, i], var)
            if i != t.nu && o != t.ny
                print("\n")
            end
        end
    end
    if iscontinuous(t)
        print(io, "\nContinuous-time transfer function model")
    else
        print(io, "\nSample Time: ")
        if t.Ts > 0
            print(io, t.Ts, " (seconds)")
        elseif t.Ts == -1
            print(io, "unspecified")
        end
        print(io, "\nDiscrete-time transfer function model")
    end
end
