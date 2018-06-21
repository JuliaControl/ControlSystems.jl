abstract type SisoTf end
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
    inputnames::Vector{String}
    outputnames::Vector{String}
    function TransferFunction{T}(matrix::Matrix{T}, Ts::Float64,
            inputnames::Vector{String}, outputnames::Vector{String}) where T<:SisoTf
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



tzero(sys::SisoTf) = error("tzero is not implemented for type $(typeof(sys))")
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
#.*(a::SisoTf, b::SisoTf) = .*(promote(a,b)...)
#.+(a::SisoTf, b::SisoTf) = .+(promote(a,b)...)
#.-(a::SisoTf, b::SisoTf) = .+(promote(a,b)...)

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
function tf{T1<:Vector, T2<:Vector}(num::VecOrMat{T1}, den::VecOrMat{T2}, Ts::Real=0; kwargs...)
    # Validate input and output dimensions match
    ny, nu = size(num, 1, 2)
    T = promote_type(T1,T2)
    if (ny, nu) != size(den, 1, 2)
        error("num and den dimensions must match")
    end
    matrix = Array{SisoRational{T}}(ny, nu)
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
function zpk{T1<:Vector,T2<:Vector}(z::VecOrMat{T1}, p::VecOrMat{T2}, k::VecOrMat, Ts::Real=0; kwargs...)
    # Validate input and output dimensions match
    ny, nu = size(z, 1, 2)
    if (ny, nu) != size(p, 1, 2) || (ny, nu) != size(k, 1, 2)
        error("s, p, and k kdimensions must match")
    end
    T = promote_type(T1,T2,Vector{Complex128})
    matrix = Array{SisoZpk{T}}(ny, nu)
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

function zpk(t1::TransferFunction)
    oldmat = t1.matrix
    vectype = primitivetype(t1)
    numbertype = promote_type(primitivetype(t1), Complex128)
    T = SisoZpk{Vector{numbertype}} # TODO: should ideally not hard code Vector
    matrix = Matrix{T}(t1.ny, t1.nu)
    for i in eachindex(oldmat)
        matrix[i] = convert(T, oldmat[i])
    end
    return TransferFunction(matrix, t1.Ts, copy(t1.inputnames), copy(t1.outputnames))
end

function tf(t1::TransferFunction)
    oldmat = t1.matrix
    T = SisoRational{Vector{primitivetype(t1)}}
    matrix = Array{T}(t1.ny, t1.nu)
    for i in eachindex(oldmat)
        matrix[i] = convert(T, oldmat[i])
    end
    return TransferFunction(matrix, t1.Ts, copy(t1.inputnames), copy(t1.outputnames))
end

@doc """ `sys = tfg(tf::LTISystem), `tfg(s::AbstractString)`, `tfg(exp::Expr)`, `tfg(::Array)`

Create generalized transfer function represented by an expression. The variable has to be `s`.

Example: `tfg("1/exp(-sqrt(s))")`, `tfg(["1/exp(-sqrt(s))"), "1/(s+1)])`, `tfg(:(s+1))`

Other uses:

`tfg(sys)`: Convert `sys` to `tfg` form.
""" ->
function tfg(t1::TransferFunction)
    oldmat = t1.matrix
    ST = SisoGeneralized{Vector{primitivetype(t1.matrix)}}
    matrix = Matrix{ST}(t1.ny, t1.nu)
    for i in eachindex(oldmat)
        matrix[i] = convert(ST, oldmat[i])
    end
    return TransferFunction(matrix, t1.Ts, copy(t1.inputnames), copy(t1.outputnames))
end

zpk(sys::TransferFunction{SisoGeneralized}) = zpk(tf(sys))

function tf(num::Vector{T1}, den::Vector{T2}, Ts::Real=0; kwargs...) where {T1,T2}
    T = promote_type(T1,T2)
    tf(reshape(Vector{T}[T.(num)], 1, 1), reshape(Vector{T}[T.(den)], 1, 1), Ts; kwargs...)
end

tf(num::Real, den::Vector, Ts::Real=0; kwargs...) = tf([num], den, Ts; kwargs...)

function zpk(z::Vector{T1}, p::Vector{T2}, k::Real, Ts::Real=0; kwargs...) where {T1,T2}
    T = if T1 == Any
        if T2 == Any
            promote_type(typeof(k), Complex128)
        else
            T2
        end
    elseif T2 == Any
        T1
    else
        promote_type(T1,T2)
    end
    zpk(reshape(Vector{T}[T.(z)], 1, 1), reshape(Vector{T}[T.(p)], 1, 1), reshape([k],1,1), Ts; kwargs...)
end

# Function for creation of static gain
function tf(gain::Array, Ts::Real=0; kwargs...)
    ny, nu = size(gain, 1, 2)
    matrix = Matrix{SisoRational{Vector{eltype(gain)}}}(ny, nu)
    for i in eachindex(gain)
        matrix[i] = SisoRational([gain[i]], [one(eltype(gain))])
    end
    kvs = Dict(kwargs)
    inputnames = validate_names(kvs, :inputnames, nu)
    outputnames = validate_names(kvs, :outputnames, ny)
    return TransferFunction(matrix, Float64(Ts), inputnames, outputnames)
end

function zpk(gain::Array, Ts::Real=0; kwargs...)
    ny, nu = size(gain, 1, 2)
    T = promote_type(Complex128,primitivetype(gain))
    matrix = Array{SisoZpk{Vector{T}}}(ny, nu)
    for o=1:ny
        for i=1:nu
            matrix[o, i] = SisoZpk(T[],T[], gain[o, i])
        end
    end
    kvs = Dict(kwargs)
    inputnames = validate_names(kvs, :inputnames, nu)
    outputnames = validate_names(kvs, :outputnames, ny)
    return TransferFunction(matrix, Float64(Ts), inputnames, outputnames)
end

tf(gain::Real, Ts::Real=0; kwargs...) = tf([gain], Ts; kwargs...)
zpk(k::Real, Ts::Real=0; kwargs...) = zpk(eltype(k)[], eltype(k)[], k, Ts; kwargs...)

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
    matrix = Matrix{SisoGeneralized}(ny, nu)
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
    rows, cols = ControlSystems.index2range(inds...)
    T = eltype(t.matrix)
    mat = Array{T}(length(rows), length(cols))
    mat[:, :] = t.matrix[rows, cols]
    innames = String[t.inputnames[i] for i in cols]
    outnames = String[t.outputnames[i] for i in rows]
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
        inputnames = fill(String(""),t1.ny)
    end
    if all(t1.outputnames .== "")
        outputnames = t2.outputnames
    elseif all(t2.outputnames .== "") || (t1.outputnames == t2.outputnames)
        outputnames = t1.outputnames
    else
        outputnames = fill(String(""),t1.nu)
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
    T = promote_type(eltype(t1.matrix), eltype(t2.matrix))
    matrix = convert(Matrix{T}, *(promote(t1.matrix,t2.matrix)...))
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
/(t1::TransferFunction, t2::TransferFunction) = t1*(one(primitivereal(t2))/t2)

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
