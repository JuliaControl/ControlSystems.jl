struct TransferFunction{TE, S<:SisoTf{T} where T} <: LTISystem{TE}
    matrix::Matrix{S}
    timeevol::TE
    function TransferFunction{TE,S}(matrix::Matrix{S}, timeevol::TE) where {S,TE}
        # Validate size of input and output names
        ny, nu = size(matrix)
        return new{TE,S}(matrix, timeevol)
    end
end
function TransferFunction(matrix::Matrix{S}, timeevol::TE) where {TE<:TimeEvolution, T<:Number, S<:SisoTf{T}}
    TransferFunction{TE, S}(matrix, timeevol)
end

# # Constructor for Discrete time system
# function TransferFunction(matrix::Matrix{S}, Ts::Number) where {T<:Number, S<:SisoTf{T}}
#     return TransferFunction(matrix, Discrete(Ts))
# end
# # Constructor for Continuous time system
# function TransferFunction(matrix::Matrix{S}) where {T<:Number,S<:SisoTf{T}}
#     return TransferFunction(matrix, Continuous())
# end

noutputs(G::TransferFunction) = size(G.matrix, 1)
ninputs(G::TransferFunction) = size(G.matrix, 2)

#####################################################################
##                          Misc. Functions                        ##
#####################################################################

## INDEXING ##
Base.ndims(::TransferFunction) = 2
Base.size(G::TransferFunction) = size(G.matrix)
Base.eltype(::Type{S}) where {S<:TransferFunction} = S
Base.zero(G::TransferFunction{TE,S}) where {TE,S} = tf(zeros(numeric_type(S), size(G)), G.timeevol) 
Base.zero(::Type{TransferFunction{TE,S}}) where {TE,S} = TransferFunction([zero(S);;], undef_sampletime(TE)) 
function Base.getproperty(G::TransferFunction, s::Symbol)
    s ∈ fieldnames(typeof(G)) && return getfield(G, s)
    if s === :ny
        return size(G, 1)
    elseif s === :nu
        return size(G, 2)
    elseif s === :Ts
        if isdiscrete(G)
            return timeevol(G).Ts
        else
            @warn "Getting time 0.0 for non-discrete systems is deprecated. Check `isdiscrete` before trying to access time."
            return 0.0
        end
        throw(ArgumentError("$(typeof(G)) has no property named $s"))
    end
end

function Base.getindex(G::TransferFunction{TE,S}, inds...) where {TE,S<:SisoTf}
    if size(inds, 1) != 2
        error("Must specify 2 indices to index TransferFunction model")
    end
    rows, cols = index2range(inds...)
    mat = Matrix{S}(
        undef,
        rows isa Colon ? noutputs(G) : length(rows),
        cols isa Colon ? ninputs(G) : length(cols)
    )
    mat[:, :] = G.matrix[rows, cols]
    return TransferFunction(mat, G.timeevol)
end

function Base.copy(G::TransferFunction)
    return TransferFunction(copy(G.matrix), G.timeevol)
end

numvec(G::TransferFunction) = map(numvec, G.matrix)
denvec(G::TransferFunction) = map(denvec, G.matrix)

numpoly(G::TransferFunction) = map(numpoly, G.matrix)
denpoly(G::TransferFunction) = map(denpoly, G.matrix)

"""
    minreal(tf::TransferFunction, eps=sqrt(eps()))

Create a minimal representation of each transfer function in `tf` by cancelling poles and zeros
will promote system to an appropriate numeric type
"""
function minreal(G::TransferFunction, eps::Real=sqrt(eps()))
    matrix = similar(G.matrix, typeof(minreal(one(first(G.matrix)), eps)))
    for i = eachindex(G.matrix)
        matrix[i] = minreal(G.matrix[i], eps)
    end
    return TransferFunction(matrix, G.timeevol)
end

"""
    isproper(tf)

Returns `true` if the `TransferFunction` is proper. This means that order(den)
>= order(num))
"""
function isproper(G::TransferFunction)
    return all(isproper(f) for f in G.matrix)
end

isrational(::TransferFunction) = true

#####################################################################
##                         Math Operators                          ##
#####################################################################

## EQUALITY ##
function ==(G1::T, G2::T) where T<:TransferFunction
    fields = (:timeevol, :ny, :nu, :matrix)
    for field in fields
        if getproperty(G1, field) != getproperty(G2, field)
            return false
        end
    end
    return true
end

==(G1::TransferFunction, G2::TransferFunction) = ==(promote(G1,G2)...)

## Approximate ##
function isapprox(G1::TransferFunction, G2::TransferFunction; kwargs...)
    G1, G2 = promote(G1, G2)
    fieldsApprox = (:timeevol, :matrix)
    for field in fieldsApprox
        if !(isapprox(getproperty(G1, field), getproperty(G2, field); kwargs...))
            return false
        end
    end
    return true
end

function isapprox(G1::Array{T}, G2::Array{S}; kwargs...) where {T<:SisoTf, S<:SisoTf}
    all(i -> isapprox(promote(G1[i], G2[i])...; kwargs...), eachindex(G1))
end

## ADDITION ##
function +(G1::TransferFunction, G2::TransferFunction)
    if size(G1) != size(G2)
        error("Systems have different shapes.")
    end
    timeevol = common_timeevol(G1,G2)
    matrix = G1.matrix + G2.matrix
    return TransferFunction(matrix, timeevol)
end

+(G::TransferFunction, n::Number) = TransferFunction(G.matrix .+ n, G.timeevol)
+(n::Number, G::TransferFunction) = +(G, n)

## SUBTRACTION ##
-(n::Number, G::TransferFunction) = TransferFunction(n .- G.matrix, G.timeevol)
-(G1::TransferFunction, G2::TransferFunction) = +(G1, -G2)
-(G::TransferFunction, n::Number) = +(G, -n)

## NEGATION ##
-(G::TransferFunction) = TransferFunction(-G.matrix, G.timeevol)

## MULTIPLICATION ##

function *(G1::TransferFunction, G2::TransferFunction)
    # Note: G1*G2 = y <- G1 <- G2 <- u
    timeevol = common_timeevol(G1,G2)
    if G1.nu != G2.ny
        if issiso(G1) || issiso(G2)
            error("G1*G2: G1 must have same number of inputs as G2 has outputs, did you intend to broadcast the multiplication?")
        else
            error("G1*G2: G1 must have same number of inputs as G2 has outputs")
        end
    end
    matrix = G1.matrix * G2.matrix
    return TransferFunction(matrix, timeevol)
end

*(G::TransferFunction, n::Number) = TransferFunction(n*G.matrix, G.timeevol)
*(n::Number, G::TransferFunction) = *(G, n)

function Base.Broadcast.broadcasted(::typeof(*), G1::TransferFunction, G2::TransferFunction)
    issiso(G1) || issiso(G2) || error("Only SISO transfer function can be broadcasted")
    # Note: G1*G2 = y <- G1 <- G2 <- u
    timeevol = common_timeevol(G1,G2)
    
    if issiso(G1) && !issiso(G2) # Check !issiso(G2) to avoid calling fill if both are siso
        G1 = append(G1 for i in 1:G2.ny)
    elseif issiso(G2)
        G2 = append(G2 for i in 1:G1.nu)
    end
    return G1 * G2
end

function Base.Broadcast.broadcasted(::typeof(*), G1::TransferFunction, M::AbstractArray{<:Number})
    issiso(G1) || error("Only SISO transfer function can be broadcasted")
    LinearAlgebra.isdiag(M) || error("Broadcasting multiplication of an LTI system with an array is only supported for diagonal arrays. If you want the system to behave like a scalar and multiply each element of the array, wrap the system in a `Ref` to indicate this, i.e., `Ref(sys) .* array`.")
    # Note: G1*G2 = y <- G1 <- G2 <- u
    timeevol = G1.timeevol
    matrix = G1.matrix .* M
    return TransferFunction(matrix, timeevol)
end

function Base.Broadcast.broadcasted(::typeof(*), M::AbstractArray{<:Number}, G1::TransferFunction)
    issiso(G1) || error("Only SISO transfer function can be broadcasted")
    LinearAlgebra.isdiag(M) || error("Broadcasting multiplication of an LTI system with an array is only supported for diagonal arrays. If you want the system to behave like a scalar and multiply each element of the array, wrap the system in a `Ref` to indicate this, i.e., `Ref(sys) .* array`.")
    # Note: G1*G2 = y <- G1 <- G2 <- u
    timeevol = G1.timeevol
    matrix = M .* G1.matrix
    return TransferFunction(matrix, timeevol)
end

function Base.Broadcast.broadcasted(::typeof(*), G1r::Base.RefValue{<:TransferFunction}, M::AbstractArray{<:Number})
    G1 = G1r[]
    issiso(G1) || error("Only SISO transfer function can be broadcasted")
    timeevol = G1.timeevol
    matrix = G1.matrix .* M
    return TransferFunction(matrix, timeevol)
end

function Base.Broadcast.broadcasted(::typeof(*), M::AbstractArray, G1r::Base.RefValue{<:TransferFunction})
    G1 = G1r[]
    issiso(G1) || error("Only SISO transfer function can be broadcasted")
    timeevol = G1.timeevol
    matrix = M .* G1.matrix
    return TransferFunction(matrix, timeevol)
end

## DIVISION ##
function /(n::Number, G::TransferFunction)
    if issiso(G)
        entry = n/G.matrix[1,1]
        matrix = fill(entry, 1, 1)
        return TransferFunction(matrix, G.timeevol)
    else
        error("MIMO TransferFunction inversion isn't implemented yet, consider converting your transfer functions to state-space form using `ss`")
    end
end
/(G::TransferFunction, n::Number) = G*(1/n)
Base.:(/)(sys1::LTISystem, sys2::TransferFunction) = *(promote(sys1, ss(1/sys2))...) # This special case is needed to properly handle improper inverse transfer function (1/s)
function /(f1::TransferFunction, f2::TransferFunction)
    if issiso(f2)
        T = numeric_type(f2)
        One = one(1/one(T)) # This gymatics is to ensure that we get floats out of integer system division. We will always get floats for ss division in the branch below, so this is required for type stability
        f1*(One/f2)
    else
        @warn "MIMO TransferFunction inversion isn't implemented yet, converting to state-space object and back. Consider converting your transfer functions to state-space form using `ss` as soon as possible."
        tf(ss(f1) / ss(f2))
    end
end



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
    elseif isdiscrete(G)
        print(io, "\nSample Time: ", G.Ts, " (seconds)")
        print(io, "\nDiscrete-time transfer function model")
    else
        print(io, "\nStatic gain transfer function model")
    end
end
