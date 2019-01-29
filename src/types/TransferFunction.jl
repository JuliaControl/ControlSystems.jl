struct TransferFunction{S<:SisoTf{T} where T} <: LTISystem
    matrix::Matrix{S}
    Ts::Float64
    nu::Int
    ny::Int
    function TransferFunction{S}(matrix::Matrix{S}, Ts::Real=0) where {S}
        # Validate size of input and output names
        ny, nu = size(matrix)
        # Validate sampling time
        if Ts < 0 && Ts != -1
            error("Ts must be either a positive number, 0
                (continuous system), or -1 (unspecified)")
        end
        return new{S}(matrix, Ts, nu, ny)
    end
end
function TransferFunction(matrix::Matrix{S}, Ts::Float64=0) where {T<:Number, S<:SisoTf{T}}
    return TransferFunction{S}(matrix, Ts)
end

noutputs(G::TransferFunction) = size(G.matrix, 1)
ninputs(G::TransferFunction) = size(G.matrix, 2)

# function TransferFunction{SisoRational{T}}(num::Matrix{Vector{T}}, den::Matrix{Vector{T}}, Ts::Real=0) where {T <: Number}
#     # Validate input and output dimensions match
#     ny, nu = size(num, 1, 2)
#     if (ny, nu) != size(den, 1, 2)
#         error("num and den dimensions must match")
#     end
#
#     matrix = Matrix{SisoRational{T}}(ny, nu) # TODO: list comprehension seems suitable here
#     for o=1:ny
#         for i=1:nu
#             matrix[o, i] = SisoRational(num[o, i], den[o, i])
#         end
#     end
#     return TransferFunction(matrix, Float64(Ts))
# end

#####################################################################
##                          Misc. Functions                        ##
#####################################################################

## INDEXING ##
Base.ndims(::TransferFunction) = 2
Base.size(G::TransferFunction) = size(G.matrix)
Base.size(G::TransferFunction, d) = size(G.matrix, d)
Base.eltype(::Type{S}) where {S<:TransferFunction} = S

function Base.getindex(G::TransferFunction{S}, inds...) where {S<:SisoTf}
    if size(inds, 1) != 2
        error("Must specify 2 indices to index TransferFunction model")
    end
    rows, cols = index2range(inds...)
    mat = Matrix{S}(undef, length(rows), length(cols))
    mat[:, :] = G.matrix[rows, cols]
    return TransferFunction(mat, G.Ts)
end

function Base.copy(G::TransferFunction)
    return TransferFunction(copy(G.matrix), G.Ts)
end

numvec(G::TransferFunction) = map(numvec, G.matrix)
denvec(G::TransferFunction) = map(denvec, G.matrix)

numpoly(G::TransferFunction) = map(numpoly, G.matrix)
denpoly(G::TransferFunction) = map(denpoly, G.matrix)

"""`tf = minreal(tf::TransferFunction, eps=sqrt(eps()))`

Create a minimial representation of each transfer function in `tf` by cancelling poles and zeros
will promote system to an appropriate numeric type"""
function minreal(G::TransferFunction, eps::Real=sqrt(eps()))
    matrix = similar(G.matrix, typeof(minreal(one(first(G.matrix)), eps)))
    for i = eachindex(G.matrix)
        matrix[i] = minreal(G.matrix[i], eps)
    end
    return TransferFunction(matrix, G.Ts)
end

"""`isproper(tf)`

Returns `true` if the `TransferFunction` is proper. This means that order(den)
>= order(num))"""
function isproper(G::TransferFunction)
    return all(isproper(f) for f in G.matrix)
end
#####################################################################
##                         Math Operators                          ##
#####################################################################

## EQUALITY ##
function ==(G1::TransferFunction, G2::TransferFunction)
    fields = [:Ts, :ny, :nu, :matrix]
    for field in fields
        if getfield(G1, field) != getfield(G2, field)
            return false
        end
    end
    return true
end

## Approximate ##
function isapprox(G1::TransferFunction, G2::TransferFunction; kwargs...)
    G1, G2 = promote(G1, G2)
    fieldsApprox = [:Ts, :matrix]
    for field in fieldsApprox
        if !(isapprox(getfield(G1, field), getfield(G2, field); kwargs...))
            return false
        end
    end
    return true
end

function isapprox(G1::Array{T}, G2::Array{S}; kwargs...) where {T<:SisoTf, S<:SisoTf}
    all(i -> isapprox(G1[i], G2[i]; kwargs...), eachindex(G1))
end

## ADDITION ##
function +(G1::TransferFunction, G2::TransferFunction)
    if size(G1) != size(G2)
        error("Systems have different shapes.")
    elseif G1.Ts != G2.Ts
        error("Sampling time mismatch")
    end

    matrix = G1.matrix + G2.matrix
    return TransferFunction(matrix, G1.Ts)
end

+(G::TransferFunction, n::Number) = TransferFunction(G.matrix .+ n, G.Ts)
+(n::Number, G::TransferFunction) = +(G, n)

## SUBTRACTION ##
-(n::Number, G::TransferFunction) = TransferFunction(n .- G.matrix, G.Ts)
-(G1::TransferFunction, G2::TransferFunction) = +(G1, -G2)
-(G::TransferFunction, n::Number) = +(G, -n)

## NEGATION ##
-(G::TransferFunction) = TransferFunction(-G.matrix, G.Ts)

## MULTIPLICATION ##
function *(G1::TransferFunction, G2::TransferFunction)
    # Note: G1*G2 = y <- G1 <- G2 <- u
    if G1.nu != G2.ny
        error("G1*G2: G1 must have same number of inputs as G2 has outputs")
    elseif G1.Ts != G2.Ts
        error("Sampling time mismatch")
    end

    matrix = G1.matrix * G2.matrix

    return TransferFunction{eltype(matrix)}(matrix, G1.Ts)
end

*(G::TransferFunction, n::Number) = TransferFunction(n*G.matrix, G.Ts)
*(n::Number, G::TransferFunction) = *(G, n)

## DIVISION ##
function /(n::Number, G::TransferFunction)
    if issiso(G)
        matrix = reshape([n/G.matrix[1,1]], 1, 1)
    else
        error("MIMO TransferFunction inversion isn't implemented yet")
    end
    return TransferFunction(matrix, G.Ts)
end
/(G::TransferFunction, n::Number) = G*(1/n)
/(G1::TransferFunction, G2::TransferFunction) = G1*(1/G2)

Base.:^(sys::TransferFunction, p::Integer) = Base.power_by_squaring(sys, p)

#####################################################################
##                        Display Functions                        ##
#####################################################################

Base.print(io::IO, G::TransferFunction) = show(io, G)

function Base.show(io::IO, G::TransferFunction)
    # Compose the name vectors
    #println(io, "TransferFunction:")
    println(io, typeof(G))
    var = iscontinuous(G) ? :s : :z
    for i=1:G.nu
        for o=1:G.ny
            if !issiso(G)
                println(io, "Input $i to output $o")
            end
                print_siso(io, G.matrix[o, i], var)
            if !(i == G.nu && o == G.ny)
                print(io, "\n")
            end
        end
    end
    if iscontinuous(G)
        print(io, "\nContinuous-time transfer function model")
    else
        print(io, "\nSample Time: ")
        if G.Ts > 0
            print(io, G.Ts, " (seconds)")
        elseif G.Ts == -1
            print(io, "unspecified")
        end
        print(io, "\nDiscrete-time transfer function model")
    end
end
