abstract SisoTf
include("polys.jl")
include("sisotf.jl")
include("sisozpk.jl")
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
    function TransferFunction{T<:SisoTf}(matrix::Matrix{T}, Ts::Float64,
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

+{T<:Real}(a::TransferFunction, b::AbstractVecOrMat{T}) = +(promote(a,b)...)

Base.promote_rule{T<:Real}(::Type{TransferFunction}, ::Type{T}) = TransferFunction
Base.promote_rule{T<:Real}(::Type{TransferFunction}, ::Union{Type{Array{T,2}},Type{Array{T,1}}}) = TransferFunction
Base.convert{T<:Real}(::Type{TransferFunction}, b::T) = tf([b], [1])

Base.promote_rule(::Type{SisoRational}, ::Type{SisoZpk}) = SisoZpk
Base.promote_rule(::Type{SisoTf}, ::Type{SisoZpk}) = SisoZpk
Base.promote_rule(::Type{SisoTf}, ::Type{SisoRational}) = SisoZpk

function Base.convert{T<:Real}(::Type{TransferFunction}, b::VecOrMat{T})
    r = Array{TransferFunction,2}(size(b,2),1)
    for j=1:size(b,2)
        r[j] = vcat(map(r->tf([r],[1]),b[:,j])...)
    end
    hcat(r...)
end

function Base.convert(::Type{SisoZpk}, sys::SisoRational)
    if length(sys.num) == 0
        return SisoZpk([],[],0)
    elseif all(sys.den == zero(sys.den))
        error("Zero denominator, this should not be possible")
    else
        return SisoZpk(roots(sys.num),roots(sys.den),sys.num[1]/sys.den[1])
    end
end

#Just default SisoTf to SisoRational
SisoTf(args...) = SisoRational(args...)
Base.convert(::Type{Control.SisoTf}, b::Real) = Base.convert(Control.SisoRational, b)
Base.zero(::Type{SisoTf}) = zero(SisoRational)
Base.zero(t::SisoTf) = zero(SisoRational)
#####################################################################
##                      SisoTf Operations                   ##
#####################################################################

#These make sure that the matrix operation below works as expected
#Base.convert(::Type{SisoTf}, b::Real) = Base.convert(SisoRational, b)
*{T<:SisoTf}(a::Array{T}, b::Real) = map(x->x*b,a)
*{T<:SisoTf}(b::Real, a::Array{T}) = map(x->x*b,a)
/{T<:SisoTf}(a::Array{T}, b::Real) = map(x->x/b,a)
+{T<:SisoTf}(a::Array{T}, b::Real) = map(x->x+b,a)
+{T<:SisoTf}(b::Real, a::Array{T}) = map(x->x+b,a)
-{T<:SisoTf}(a::Array{T}, b::Real) = map(x->x-b,a)
-{T<:SisoTf}(b::Real, a::Array{T}) = map(x->b-x,a)
-{T<:SisoTf}(a::Array{T})          = map(x-> -x,a)

#Operations with different types of Siso functions
*(a::SisoTf, b::SisoTf)  = *(promote(a,b)...)
+(a::SisoTf, b::SisoTf)  = +(promote(a,b)...)
-(a::SisoTf, b::SisoTf)  = -(promote(a,b)...)
.*(a::SisoTf, b::SisoTf) = .*(promote(a,b)...)
.+(a::SisoTf, b::SisoTf) = .+(promote(a,b)...)
.-(a::SisoTf, b::SisoTf) = .+(promote(a,b)...)
#####################################################################
##                      Constructor Functions                      ##
#####################################################################

function tf{T<:Vector}(num::VecOrMat{T}, den::VecOrMat{T}, Ts::Real=0; kwargs...)
    # Validate input and output dimensions match
    ny, nu = size(num, 1, 2)
    if (ny, nu) != size(den, 1, 2)
        error("num and den dimensions must match")
    end
    matrix = Array(SisoRational, ny, nu)
    for o=1:ny
        for i=1:nu
            matrix[o, i] = SisoRational(num[o, i], den[o, i])
        end
    end
    kvs = Dict(kwargs)
    inputnames = validate_names(kvs, :inputnames, nu)
    outputnames = validate_names(kvs, :outputnames, ny)
    return TransferFunction(matrix, Float64(Ts), inputnames, outputnames)
end

function zpk{T<:Vector}(z::VecOrMat{T}, p::VecOrMat{T}, k::VecOrMat, Ts::Real=0; kwargs...)
    # Validate input and output dimensions match
    ny, nu = size(z, 1, 2)
    if (ny, nu) != size(p, 1, 2) || (ny, nu) != size(k, 1, 2)
        error("s, p, and k kdimensions must match")
    end
    matrix = Array(SisoZpk, ny, nu)
    for o=1:ny
        for i=1:nu
            matrix[o, i] = SisoZpk(z[o, i], p[o, i], k[o, i])
        end
    end
    kvs = Dict(kwargs)
    inputnames = validate_names(kvs, :inputnames, nu)
    outputnames = validate_names(kvs, :outputnames, ny)
    return TransferFunction(matrix, Float64(Ts), inputnames, outputnames)
end

function zpk(tf::TransferFunction)
    matrix = Array(SisoZpk, tf.ny, tf.nu)
    for o=1:tf.ny
        for i=1:tf.nu
            matrix[o, i] = convert(SisoZpk, tf.matrix[o, i])
        end
    end
    tf.matrix = matrix
    return tf
end

tf(num::Vector, den::Vector, args...) =
tf(reshape(Vector[num], 1, 1), reshape(Vector[den], 1, 1), args...)

tf(num::Real, den::Vector, args...) = tf([num], den, args...)

zpk(z::Vector, p::Vector, k::Real, args...) =
zpk(reshape(Vector[z], 1, 1), reshape(Vector[p], 1, 1), reshape([k],1,1), args...)

# Function for creation of static gain
function tf(gain::Array, Ts::Real=0; kwargs...)
    ny, nu = size(gain, 1, 2)
    matrix = Array(SisoRational, ny, nu)
    for o=1:ny
        for i=1:nu
            matrix[o, i] = SisoRational([gain[o, i]], [1])
        end
    end
    kvs = Dict(kwargs)
    inputnames = validate_names(kvs, :inputnames, nu)
    outputnames = validate_names(kvs, :outputnames, ny)
    return TransferFunction(matrix, Float64(Ts), inputnames, outputnames)
end
tf(gain::Real, Ts::Real=0; kwargs...) = tf([gain], Ts, kwargs...)

zpk(gain::Array, Ts::Real=0; kwargs...) = zpk(tf(gain, Ts; kwargs...))

# Function for creation of 's' or 'z' var
function tf(var::AbstractString)
    var != "s" && error("var must be 's' for continuous time tf.")
    return tf([1, 0], [1])
end
function tf(var::AbstractString, Ts::Real)
    var != "z" && error("var must be 'z' for discrete time tf.")
    Ts == 0 && error("Ts must not be 0 for discrete time tf.")
    return tf([1, 0], [1], Ts)
end

zpk(var::AbstractString) = zpk(tf(var))
zpk(var::AbstractString, Ts::Real) = zpk(tf(var, Ts))

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
    mat = Array(typeof(t.matrix).parameters[1], length(rows), length(cols))
    mat[:, :] = t.matrix[rows, cols]
    innames = length(cols) > 1 ? collect(t.inputnames[cols]) : [t.inputnames[cols]];
    outnames = length(rows) > 1 ? collect(t.outputnames[rows]) : [t.inputnames[rows]];
    return TransferFunction(mat, t.Ts, innames, outnames)
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
    fields = [:Ts, :ny, :nu, :inputnames, :outputnames, :matrix]
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
                print_siso(io, t.matrix[o, i], var)
            if !(i == t.nu && o == t.ny)
                print(io, "\n")
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
