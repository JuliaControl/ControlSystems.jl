#####################################################################
##                      Data Type Declarations                     ##
#####################################################################

struct ExtendedStateSpace{T, MT<:AbstractMatrix{T}} <: LTISystem
    A::MT
    B1::MT
    B2::MT
    C1::MT
    C2::MT
    D11::MT
    D12::MT
    D21::MT
    D22::MT
    Ts::Float64
    nx::Int
    nw::Int
    nu::Int
    nz::Int
    ny::Int
    function ExtendedStateSpace{T, MT}(A::MT, B1::MT, B2::MT,C1::MT,C2::MT,
        D11::MT, D12::MT, D21::MT, D22::MT, Ts::Float64) where {T, MT <: AbstractMatrix{T}}
        nx = size(A, 1)
        nw = size(B1, 2)
        nu = size(B2, 2)
        nz = size(C1, 1)
        ny = size(C2, 1)

        if size(A, 2) != nx && nx != 0
            error("A must be square")
        elseif size(B1, 1) != nx
            error("B1 must have the same row size as A")
          elseif size(B2, 1) != nx
              error("B2 must have the same row size as A")
        elseif size(C1, 2) != nx
            error("C1 must have the same column size as A")
        elseif size(C2, 2) != nx
            error("C2 must have the same column size as A")
        elseif nw != size(D11, 2)
            error("D11 must have the same column size as B1")
        elseif nw != size(D21, 2)
            error("D12 must have the same column size as B1")
        elseif nu != size(D12, 2)
            error("D12 must have the same column size as B2")
        elseif nu != size(D22, 2)
            error("D22 must have the same column size as B2")
        elseif nz != size(D11, 1)
            error("D11 must have the same row size as C1")
        elseif nz != size(D12, 1)
            error("D12 must have the same row size as C1")
        elseif ny != size(D21, 1)
            error("D11 must have the same row size as C12")
        elseif ny != size(D22, 1)
            error("D12 must have the same row size as C2")
        end

        # Validate sampling time
        if Ts < 0 && Ts != -1
            error("Ts must be either a positive number, 0
                   (continuous system), or -1 (unspecified)")
        end
        new{T, MT}(A, B1, B2, C1, C2, D11, D12, D21, D22, Ts, nx, nw, nu, nz, ny)
    end
end

function ExtendedStateSpace{T,MT}(A::AbstractArray, B1::AbstractArray,  B2::AbstractArray, C1::AbstractArray, C2::AbstractArray, D11::AbstractArray, D12::AbstractArray, D21::AbstractArray, D22::AbstractArray, Ts::Real) where {T, MT <: AbstractMatrix{T}}
        return ExtendedStateSpace{T,Matrix{T}}(MT(to_matrix(T, A)), MT(to_matrix(T, B1)),
            MT(to_matrix(T, B2)), MT(to_matrix(T, C1)), MT(to_matrix(T, C2)),
            MT(to_matrix(T, D11)), MT(to_matrix(T, D12)), MT(to_matrix(T, D21)),
            MT(to_matrix(T, D22)), Float64(Ts))
end

function ExtendedStateSpace(A::AbstractArray, B1::AbstractArray, B2::AbstractArray, C1::AbstractArray, C2::AbstractArray, D11::AbstractArray, D12::AbstractArray, D21::AbstractArray, D22::AbstractArray, Ts::Real)
        # TODO: change back in 0.7 T = promote_type(eltype(A),eltype(B),eltype(C),eltype(D))
        TBC = promote_type(promote_type(eltype(B1),eltype(B2)),promote_type(eltype(C1),eltype(C2)))
        TD  = promote_type(promote_type(eltype(D11),eltype(D12)),promote_type(eltype(D21),eltype(D22)))
        T   = promote_type(promote_type(TBC,TD),eltype(A))
        @assert (typeof(to_matrix(T, A)) == typeof(to_matrix(T, B1)) == typeof(to_matrix(T, B2)) == typeof(to_matrix(T, C1)) == typeof(to_matrix(T, C2)) == typeof(to_matrix(T, D11)) == typeof(to_matrix(T, D12)) == typeof(to_matrix(T, D21)))
        return ExtendedStateSpace{T,Matrix{T}}(to_matrix(T, A), to_matrix(T, B1),
            to_matrix(T, B2), to_matrix(T, C1), to_matrix(T, C2), to_matrix(T, D11),
            to_matrix(T, D12), to_matrix(T, D21), to_matrix(T, D22), Float64(Ts))
end

# Getter functions
get_A(sys::ExtendedStateSpace)   = sys.A
get_B1(sys::ExtendedStateSpace)  = sys.B1
get_B2(sys::ExtendedStateSpace)  = sys.B2
get_B(sys::ExtendedStateSpace)   = [sys.B1 sysB2]
get_C1(sys::ExtendedStateSpace)  = sys.C1
get_C2(sys::ExtendedStateSpace)  = sys.C2
get_C(sys::ExtendedStateSpace)   = [sys.C1; sys.C2]
get_D11(sys::ExtendedStateSpace) = sys.D11
get_D12(sys::ExtendedStateSpace) = sys.D12
get_D21(sys::ExtendedStateSpace) = sys.D21
get_D22(sys::ExtendedStateSpace) = sys.D22
get_D(sys::ExtendedStateSpace)   = [sys.D11 sys.D12 ; sys.D21 sys.D22]

get_Ts(sys::ExtendedStateSpace) = sys.Ts

ssdata(sys::ExtendedStateSpace) = get_A(sys), get_B1(sys), get_B2(sys), get_C1(sys), get_C2(sys), get_D11(sys), get_D12(sys), get_D21(sys), get_D22(sys)

# Funtions for number of intputs, outputs and states
ninputs(sys::ExtendedStateSpace) = size(get_D(sys), 2)
noutputs(sys::ExtendedStateSpace) = size(get_D(sys), 1)
nstates(sys::ExtendedStateSpace) = size(get_A(sys), 1)

#####################################################################
##                         Math Operators                          ##
#####################################################################

## EQUALITY ##
function ==(sys1::ExtendedStateSpace, sys2::ExtendedStateSpace)
    return all(getfield(sys1, f) == getfield(sys2, f) for f in fieldnames(ExtendedStateSpace))
end

## Approximate ##
function isapprox(sys1::ExtendedStateSpace, sys2::ExtendedStateSpace)
    return all(getfield(sys1, f) ≈ getfield(sys2, f) for f in fieldnames(ExtendedStateSpace))
end

## ADDITION ##
# not sure how to best handle this yet

## SUBTRACTION ##
# not sure how to best handle this yet

## NEGATION ##
# not sure how to best handle this yet

## MULTIPLICATION ##
# not sure how to best handle this yet

## DIVISION ##
# not sure how to best handle this yet

#####################################################################
##                       Indexing Functions                        ##
#####################################################################
Base.ndims(::ExtendedStateSpace) = 2 # NOTE: Also for SISO systems?
Base.size(sys::ExtendedStateSpace) = (noutputs(sys), ninputs(sys)) # NOTE: or just size(get_D(sys))
Base.size(sys::ExtendedStateSpace, d) = d <= 2 ? size(sys)[d] : 1
Base.eltype(::Type{S}) where {S<:ExtendedStateSpace} = S

function Base.getindex(sys::ExtendedStateSpace, inds...)
    if size(inds, 1) != 2
        error("Must specify 2 indices to index statespace model")
    end
    rows, cols = index2range(inds...) # FIXME: ControlSystems.index2range(inds...)
    return ExtendedStateSpace(copy(sys.A), sys.B[:, cols], sys.C[rows, :], sys.D[rows, cols], sys.Ts)
end

#####################################################################
##                        Display Functions                        ##
#####################################################################

Base.print(io::IO, sys::ExtendedStateSpace) = show(io, sys)

function Base.show(io::IO, sys::ExtendedStateSpace)
    # Compose the name vectors
    println(io, typeof(sys))
    if nstates(sys) > 0
        println(io, "A = \n", _string_mat_with_headers(sys.A))
        println(io, "B1 = \n", _string_mat_with_headers(sys.B1))
        println(io, "B2 = \n", _string_mat_with_headers(sys.B2))
        println(io, "C1 = \n", _string_mat_with_headers(sys.C1))
        println(io, "C2 = \n", _string_mat_with_headers(sys.C2))
        println(io, "D11 = \n", _string_mat_with_headers(sys.D11))
        println(io, "D12 = \n", _string_mat_with_headers(sys.D12))
        println(io, "D21 = \n", _string_mat_with_headers(sys.D21))
        println(io, "D22 = \n", _string_mat_with_headers(sys.D22))
    else
        println(io, "The extended statespece model has no states..!")
    end

    # Print sample time
    if sys.Ts > 0
        println(io, "Sample Time: ", sys.Ts, " (seconds)")
    elseif sys.Ts == -1
        println(io, "Sample Time: unspecified")
    end

    # Print model type
    if nstates(sys) == 0
        print(io, "Static gain")
    elseif iscontinuous(sys)
        print(io, "Continuous-time extended state-space model")
    else
        print(io, "Discrete-time extended state-space model")
    end
end
