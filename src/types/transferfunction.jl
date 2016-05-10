abstract SisoTf
include("polys.jl")
include("sisotf.jl")
include("sisozpk.jl")
include("sisogeneralized.jl")
#####################################################################
##                      Data Type Declarations                     ##
#####################################################################

type TransferFunction{S<:SisoTf} <: LTISystem
    matrix::Matrix{S}
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
        return new{T}(matrix, Ts, nu, ny, inputnames, outputnames)
    end
end
TransferFunction{T<:SisoTf}(matrix::Matrix{T}, args...) = TransferFunction{T}(matrix, args...)

+{T<:Real}(a::TransferFunction, b::AbstractVecOrMat{T}) = +(promote(a,b)...)

Base.promote_rule{S<:SisoTf,T<:Real}(::Type{TransferFunction{S}}, ::Type{T}) = TransferFunction{S}
Base.promote_rule{S<:SisoTf,T<:Real}(::Type{TransferFunction{S}}, ::Union{Type{Array{T,2}},Type{Array{T,1}}}) = TransferFunction{S}

Base.convert{T<:Real}(::Type{TransferFunction}, b::T) = tf([b])
Base.convert{T<:Real}(::Type{TransferFunction{SisoRational}}, b::T) = tf(b)
Base.convert{T<:Real}(::Type{TransferFunction{SisoZpk}}, b::T) = zpk(b)
Base.convert{T<:Real}(::Type{TransferFunction{SisoGeneralized}}, b::T) = tfg(b)

Base.convert(::Type{TransferFunction{SisoZpk}}, s::TransferFunction) = zpk(s)
Base.convert(::Type{TransferFunction{SisoRational}}, s::TransferFunction) = tf(s)
Base.convert(::Type{TransferFunction{SisoGeneralized}}, s::TransferFunction) = tfg(s)

Base.promote_rule(::Type{TransferFunction{SisoRational}}, ::Type{TransferFunction{SisoZpk}}) = TransferFunction{SisoZpk}
Base.promote_rule{T<:SisoTf}(::Type{TransferFunction{T}}, ::Type{TransferFunction{SisoGeneralized}}) = TransferFunction{SisoGeneralized}
Base.promote_rule(::Type{SisoRational}, ::Type{SisoZpk}) = SisoZpk
Base.promote_rule{T<:SisoTf}(::Type{T}, ::Type{SisoGeneralized}) = SisoGeneralized

function Base.convert{T<:Real}(::Type{TransferFunction}, b::VecOrMat{T})
    r = Array{TransferFunction,2}(size(b,2),1)
    for j=1:size(b,2)
        r[j] = vcat(map(k->convert(TransferFunction,k),b[:,j])...)
    end
    hcat(r...)
end

function Base.convert(::Type{SisoZpk}, sys::SisoRational)
    if length(numvec(sys)) == 0
        return SisoZpk([],[],0)
    elseif all(denvec(sys) == zero(denvec(sys)))
        error("Zero denominator, this should not be possible")
    else
        return SisoZpk(roots(numpoly(sys)),roots(denpoly(sys)),numvec(sys)[1]/denvec(sys)[1])
    end
end

function Base.convert(::Type{SisoRational}, sys::SisoZpk)
    return SisoRational(numpoly(sys), denpoly(sys))
end

Base.convert(::Type{SisoGeneralized}, sys::SisoRational) = SisoGeneralized(sprint(print_compact, sys))
Base.convert(::Type{SisoGeneralized}, sys::SisoZpk) = convert(SisoGeneralized, convert(SisoRational, sys))

Base.convert(::Type{SisoRational}, sys::SisoGeneralized) = SisoRational(sys.expr)
Base.convert(::Type{SisoZpk}, sys::SisoGeneralized) = convert(SisoZpk, SisoRational(sys.expr))
Base.convert(::Type{TransferFunction{SisoZpk}}, s::TransferFunction) = zpk(s)

#Just default SisoTf to SisoRational
SisoTf(args...) = SisoRational(args...)
Base.convert(::Type{ControlSystems.SisoTf}, b::Real) = Base.convert(ControlSystems.SisoRational, b)
Base.zero(::Type{SisoTf}) = zero(SisoRational)
Base.zero(::SisoTf) = zero(SisoRational)

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

==(a::SisoTf, b::SisoTf) = ==(promote(a,b)...)
!=(a::SisoTf, b::SisoTf) = !(a==b)
isapprox(a::SisoTf, b::SisoTf; kwargs...) = isapprox(promote(a,b)...; kwargs...)
#####################################################################
##                      Constructor Functions                      ##
#####################################################################

@doc """ `sys = tf(num, den, Ts=0; kwargs...), sys = tf(gain, Ts=0; kwargs...)`

Create transfer function as a fraction of polynomials:

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

@doc """ `zpk(gain, Ts=0; kwargs...), zpk(num, den, k, Ts=0; kwargs...), zpk(sys)`

Create transfer function on zero pole gain form. The numerator and denominator are represented by their poles and zeros.

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
    oldmat = tf.matrix
    matrix = Array(SisoZpk, tf.ny, tf.nu)
    for i in eachindex(oldmat)
        matrix[i] = convert(SisoZpk, oldmat[i])
    end
    return TransferFunction(matrix, tf.Ts, copy(tf.inputnames), copy(tf.outputnames))
end

function tf(tf::TransferFunction)
    oldmat = tf.matrix
    matrix = Array(SisoRational, tf.ny, tf.nu)
    for i in eachindex(oldmat)
        matrix[i] = convert(SisoRational, oldmat[i])
    end
    return TransferFunction(matrix, tf.Ts, copy(tf.inputnames), copy(tf.outputnames))
end

@doc """ `sys = tfg(tf::LTISystem), `tfg(s::AbstractString)`, `tfg(exp::Expr)`, `tfg(::Array)`

Create generalized transfer function represented by an expression. The variable has to be `s`.

Example: `tfg("1/exp(-sqrt(s))")`, `tfg(["1/exp(-sqrt(s))"), "1/(s+1)])`, `tfg(:(s+1))`

Other uses:

`tfg(sys)`: Convert `sys` to `tfg` form.
""" ->
function tfg(tf::TransferFunction)
    oldmat = tf.matrix
    matrix = Array(SisoGeneralized, tf.ny, tf.nu)
    for i in eachindex(oldmat)
        matrix[i] = convert(SisoGeneralized, oldmat[i])
    end
    return TransferFunction(matrix, tf.Ts, copy(tf.inputnames), copy(tf.outputnames))
end

zpk(sys::TransferFunction{SisoGeneralized}) = zpk(tf(sys))

tf(num::Vector, den::Vector, Ts::Real=0; kwargs...) =
    tf(reshape(Vector[num], 1, 1), reshape(Vector[den], 1, 1), Ts; kwargs...)

tf(num::Real, den::Vector, Ts::Real=0; kwargs...) = tf([num], den, Ts; kwargs...)

zpk(z::Vector, p::Vector, k::Real, Ts::Real=0; kwargs...) =
    zpk(reshape(Vector[z], 1, 1), reshape(Vector[p], 1, 1), reshape([k],1,1), Ts; kwargs...)

# Function for creation of static gain
function tf(gain::Array, Ts::Real=0; kwargs...)
    ny, nu = size(gain, 1, 2)
    matrix = Array(SisoRational, ny, nu)
    for i in eachindex(gain)
        matrix[i] = SisoRational([gain[i]], [1])
    end
    kvs = Dict(kwargs)
    inputnames = validate_names(kvs, :inputnames, nu)
    outputnames = validate_names(kvs, :outputnames, ny)
    return TransferFunction(matrix, Float64(Ts), inputnames, outputnames)
end

function zpk(gain::Array, Ts::Real=0; kwargs...)
    ny, nu = size(gain, 1, 2)
    matrix = Array(SisoZpk, ny, nu)
    for o=1:ny
        for i=1:nu
            matrix[o, i] = SisoZpk([],[], gain[o, i])
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

function tfg(systems::Array, Ts::Real=0; kwargs...)
    ny, nu = size(systems, 1, 2)
    matrix = Array(SisoGeneralized, ny, nu)
    for o=1:ny
        for i=1:nu
            matrix[o, i] = SisoGeneralized(systems[o, i])
        end
    end
    kvs = Dict(kwargs)
    inputnames = validate_names(kvs, :inputnames, nu)
    outputnames = validate_names(kvs, :outputnames, ny)
    return TransferFunction(matrix, Float64(Ts), inputnames, outputnames)
end

tfg(var::Union{AbstractString,ExprLike}, Ts=0; kwargs...) = tfg([var], Ts; kwargs...)
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

@doc """`tf = minreal(tf::TransferFunction, eps=sqrt(eps()))`

Create a minimial representation of each transfer function in `tf` by cancelling poles and zeros """ ->
function minreal(t::TransferFunction, eps::Real=sqrt(eps()))
    matrix = similar(t.matrix)
    for i = eachindex(t.matrix)
        matrix[i] = minreal(t.matrix[i], eps)
    end
    return TransferFunction(matrix, t.Ts, copy(t.inputnames), copy(t.outputnames))
end

# Fallback for new Siso Types
tzero{T<:SisoTf}(::T) = error("tzero is not implemented for type $T")
pole{T<:SisoTf}(::T) = error("pole is not implemented for type $T")
numvec{T<:SisoTf}(::T) = error("numvec is not implemented for type $T")
denvec{T<:SisoTf}(::T) = error("denvec is not implemented for type $T")
numpoly{T<:SisoTf}(::T) = error("numpoly is not implemented for type $T")
denpoly{T<:SisoTf}(::T) = error("denpoly is not implemented for type $T")

"""
`numpolys = numpoly(tf::TransferFunction)`

Get an `Array` of size `(ny,nu)` containing `Poly`s representing the
numerators of `tf` from each input to output.

The numerators from `numpoly` divided by the denominators in `denpoly`
contains a full respresentation of `tf`.

See also `denpoly`, `numvec`, `denvec`.
"""
numpoly(t::TransferFunction) = map(numpoly, t.matrix)

"""
`denpolys = denpoly(tf::TransferFunction)`

Get an `Array` of size `(ny,nu)` containing `Poly`s representing the
denominators of `tf` from each input to output.

The numerators from `numpoly` divided by the denominators in `denpoly`
contains a full respresentation of `tf`.

See also `numpoly`, `numvec`, `denvec`.
"""
denpoly(t::TransferFunction) = map(denpoly, t.matrix)

"""
`numvecs = numvec(tf::TransferFunction)`

Get an `Array` of size `(ny,nu)` containing `Vector`s representing the
numerators of `tf` from each input to output.

The numerators from `numpoly` divided by the denominators in `denpoly`
contains a full respresentation of `tf`.

The polynomials are represented by the coefficients, from largest to smallest
exponent, ending with the constant term.

See also `numpoly`, `denpoly`, `denvec`.
"""
numvec(t::TransferFunction) = map(numvec, t.matrix)

"""
`denvecs = denvec(tf::TransferFunction)`

Get an `Array` of size `(ny,nu)` containing `Vector`s representing the
denominators of `tf` from each input to output.

The numerators from `numpoly` divided by the denominators in `denpoly`
contains a full respresentation of `tf`.

The polynomials are represented by the coefficients, from largest to smallest
exponent, ending with the constant term.

See also `numpoly`, `denpoly`, `numvec`.
"""
denvec(t::TransferFunction) = map(denvec, t.matrix)
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

## Approximate ##
function isapprox(t1::TransferFunction, t2::TransferFunction; kwargs...)
    t1, t2 = promote(t1, t2)
    fieldsApprox = [:Ts, :ny, :nu, :matrix]
    fieldsEqual = [:inputnames, :outputnames]
    for field in fieldsApprox
        if !(isapprox(getfield(t1, field), getfield(t2, field); kwargs...))
            return false
        end
    end
    for field in fieldsEqual
        if getfield(t1, field) != getfield(t2, field)
            return false
        end
    end
    return true
end

function isapprox{T<:SisoTf, S<:SisoTf}(t1::Array{T}, t2::Array{S}; kwargs...)
    reduce(&, [isapprox(t1[i], t2[i]; kwargs...) for i in eachindex(t1)])
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
