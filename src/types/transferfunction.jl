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
Base.promote_rule(::Type{SisoTf}, ::Type{SisoRational}) = SisoRational

function Base.convert{T<:Real}(::Type{TransferFunction}, b::VecOrMat{T})
    r = Array{TransferFunction,2}(size(b,2),1)
    for j=1:size(b,2)
        r[j] = vcat(map(s->tf([s],[1]),b[:,j])...)
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

function Base.convert(::Type{SisoRational}, sys::SisoZpk)
    num = prod(zp2polys(sys.z))*sys.k
    den = prod(zp2polys(sys.p))
    return SisoRational(num, den)
end

#Just default SisoTf to SisoRational
SisoTf(args...) = SisoRational(args...)
Base.convert(::Type{ControlSystems.SisoTf}, b::Real) = Base.convert(ControlSystems.SisoRational, b)
Base.zero(::Type{SisoTf}) = zero(SisoRational)
Base.zero(::SisoTf) = zero(SisoRational)

tzero(sys::SisoTf) = roots(sys.num)
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

@doc """ `tf(num, den, Ts=0; kwargs...), tf(gain, Ts=0; kwargs...)` Create transfer function as a fraction of polynomials:

`sys = numerator/denominator`

`num`: the coefficients of the numerator polynomial. Either scalar or vector to create SISO systems
or an array of vectors to create MIMO system.

`den`: the coefficients of the denominator polynomial. Either vector to create SISO systems
or an array of vectors to create MIMO system.

`Ts`: Sample time or `0` for continuous system.

`kwargs`: `inputnames`, `outputnames`: Arrays of strings representing the inputs and outputs.

Other uses:

`tf(sys)`: Convert `sys` to `tf` form.

`tf("s")`, `tf("z")`: Create the continous transferfunction `s`.""" ->
function tf{T<:Vector, S<:Vector}(num::VecOrMat{T}, den::VecOrMat{S}, Ts::Real=0; kwargs...)
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

@doc """ `zpk(gain, Ts=0; kwargs...), zpk(num, den, k, Ts=0; kwargs...), zpk(sys)` Create transfer
function on zero pole gain form. The numerator and denominator are represented by their poles and zeros.

`sys = k*numerator/denominator`

`num`: the roots of the numerator polynomial. Either scalar or vector to create SISO systems
or an array of vectors to create MIMO system.

`den`: the roots of the denominator polynomial. Either vector to create SISO systems
or an array of vectors to create MIMO system.

`k`: The gain of the system. Obs, this is not the same as `dcgain`.

`Ts`: Sample time or `0` for continuous system.

`kwargs`: `inputnames`, `outputnames`: Arrays of strings representing the inputs and outputs.

Other uses:

`tf(sys)`: Convert `sys` to `tf` form.

`tf("s")`: Create the transferfunction `s`.""" ->
function zpk{T<:Vector,S<:Vector}(z::VecOrMat{T}, p::VecOrMat{S}, k::VecOrMat, Ts::Real=0; kwargs...)
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
    tf = copy(tf)
    matrix = Array(SisoZpk, tf.ny, tf.nu)
    for o=1:tf.ny
        for i=1:tf.nu
            matrix[o, i] = convert(SisoZpk, tf.matrix[o, i])
        end
    end
    tf.matrix = matrix
    return tf
end

function tf(s::TransferFunction)
    s = copy(s)
    matrix = Array(SisoRational, s.ny, s.nu)
    for o=1:s.ny
        for i=1:s.nu
            matrix[o, i] = convert(SisoRational, s.matrix[o, i])
        end
    end
    s.matrix = matrix
    return s
end

tf(num::Vector, den::Vector, Ts::Real=0; kwargs...) =
    tf(reshape(Vector[num], 1, 1), reshape(Vector[den], 1, 1), Ts; kwargs...)

tf(num::Real, den::Vector, Ts::Real=0; kwargs...) = tf([num], den, Ts; kwargs...)

zpk(z::Vector, p::Vector, k::Real, Ts::Real=0; kwargs...) =
    zpk(reshape(Vector[z], 1, 1), reshape(Vector[p], 1, 1), reshape([k],1,1), Ts; kwargs...)

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
tf(gain::Real, Ts::Real=0; kwargs...) = tf([gain], Ts; kwargs...)
zpk(k::Real, Ts::Real=0; kwargs...) = zpk([], [], k, Ts; kwargs...)

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
    mat = Array(eltype(t.matrix), length(rows), length(cols))
    mat[:, :] = t.matrix[rows, cols]
    innames = UTF8String[t.inputnames[i] for i in cols]
    outnames = UTF8String[t.outputnames[i] for i in rows]
    return TransferFunction(mat, t.Ts, innames, outnames)
end

function Base.copy(t::TransferFunction)
    matrix = copy(t.matrix)
    inputnames = copy(t.inputnames)
    outputnames = copy(t.outputnames)
    return TransferFunction(matrix, t.Ts, inputnames, outputnames)
end

function minreal(t::TransferFunction, eps::Real=sqrt(eps()))
    matrix = similar(t.matrix)
    for o=1:t.ny
        for i=1:t.nu
            matrix[o, i] = minreal(t.matrix[o, i], eps)
        end
    end
    return TransferFunction(matrix, t.Ts, copy(t.inputnames), copy(t.outputnames))
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
        inputnames = fill(UTF8String(""),t1.ny)
    end
    if all(t1.outputnames .== "")
        outputnames = t2.outputnames
    elseif all(t2.outputnames .== "") || (t1.outputnames == t2.outputnames)
        outputnames = t1.outputnames
    else
        outputnames = fill(UTF8String(""),t1.nu)
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
    tftype = iscontinuous(t) ? :s : :z
    for i=1:t.nu
        for o=1:t.ny
            if !issiso(t)
                println(io, inputs[i], " to ", outputs[o])
            end
                print_siso(io, t.matrix[o, i], tftype)
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
