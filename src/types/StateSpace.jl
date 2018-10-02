#####################################################################
##                      Data Type Declarations                     ##
#####################################################################

struct StateSpace{T, MT<:AbstractMatrix{T}} <: LTISystem
    A::MT
    B::MT
    C::MT
    D::MT
    Ts::Float64
    nx::Int
    nu::Int
    ny::Int
    function StateSpace{T, MT}(A::MT, B::MT,
            C::MT, D::MT, Ts::Float64) where {T, MT <: AbstractMatrix{T}}
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

        # Validate sampling time
        if Ts < 0 && Ts != -1
            error("Ts must be either a positive number, 0
                   (continuous system), or -1 (unspecified)")
        end
        new{T, MT}(A, B, C, D, Ts, nx, nu, ny)
    end
end

function StateSpace{T,MT}(A::AbstractArray, B::AbstractArray, C::AbstractArray, D::AbstractArray, Ts::Real) where {T, MT <: AbstractMatrix{T}}
        return StateSpace{T,Matrix{T}}(MT(to_matrix(T, A)), MT(to_matrix(T, B)), MT(to_matrix(T, C)),
            MT(to_matrix(T, D)), Float64(Ts))
end

function StateSpace(A::AbstractArray, B::AbstractArray, C::AbstractArray, D::AbstractArray, Ts::Real)
        # TODO: change back in 0.7 T = promote_type(eltype(A),eltype(B),eltype(C),eltype(D))
        T = promote_type(promote_type(eltype(A),eltype(B)), promote_type(eltype(C),eltype(D)))
        @assert (typeof(to_matrix(T, A)) == typeof(to_matrix(T, B)) == typeof(to_matrix(T, C)) == typeof(to_matrix(T, D)))
        return StateSpace{T,Matrix{T}}(to_matrix(T, A), to_matrix(T, B), to_matrix(T, C),
            to_matrix(T, D), Float64(Ts))
end

# Getter functions
get_A(sys::StateSpace) = sys.A
get_B(sys::StateSpace) = sys.B
get_C(sys::StateSpace) = sys.C
get_D(sys::StateSpace) = sys.D

get_Ts(sys::StateSpace) = sys.Ts

ssdata(sys::StateSpace) = get_A(sys), get_B(sys), get_C(sys), get_D(sys)

# Funtions for number of intputs, outputs and states
ninputs(sys::StateSpace) = size(get_D(sys), 2)
noutputs(sys::StateSpace) = size(get_D(sys), 1)
nstates(sys::StateSpace) = size(get_A(sys), 1)

#####################################################################
##                         Math Operators                          ##
#####################################################################

## EQUALITY ##
function ==(sys1::StateSpace, sys2::StateSpace)
    return all(getfield(sys1, f) == getfield(sys2, f) for f in fieldnames(StateSpace))
end

## Approximate ##
function isapprox(sys1::StateSpace, sys2::StateSpace)
    return all(getfield(sys1, f) ≈ getfield(sys2, f) for f in fieldnames(StateSpace))
end

## ADDITION ##
function +(s1::StateSpace{T,MT}, s2::StateSpace{T,MT}) where {T, MT}
    #Ensure systems have same dimensions
    if size(s1) != size(s2)
        error("Systems have different shapes.")
    elseif s1.Ts != s2.Ts
        error("Sampling time mismatch")
    end

    A = [s1.A                   fill(zero(T), nstates(s1), nstates(s2));
         fill(zero(T), nstates(s2), nstates(s1))        s2.A]
    B = [s1.B ; s2.B]
    C = [s1.C s2.C;]
    D = [s1.D + s2.D;]

    return StateSpace(A, B, C, D, s1.Ts)
end

+(sys::StateSpace, n::Number) = StateSpace(sys.A, sys.B, sys.C, sys.D .+ n, sys.Ts)
+(n::Number, sys::StateSpace) = +(sys, n)

## SUBTRACTION ##
-(sys1::StateSpace, sys2::StateSpace) = +(sys1, -sys2)
-(sys::StateSpace, n::Number) = +(sys, -n)
-(n::Number, sys::StateSpace) = +(-sys, n)

## NEGATION ##
-(sys::StateSpace) = StateSpace(sys.A, sys.B, -sys.C, -sys.D, sys.Ts)

## MULTIPLICATION ##
function *(sys1::StateSpace{T,MT}, sys2::StateSpace{T,MT}) where {T, MT}
    #Check dimension alignment
    #Note: sys1*sys2 = y <- sys1 <- sys2 <- u
    if sys1.nu != sys2.ny
        error("sys1*sys2: sys1 must have same number of inputs as sys2 has outputs")
    elseif sys1.Ts != sys2.Ts
        error("Sampling time mismatch")
    end

    A = [sys1.A    sys1.B*sys2.C;
         fill(zero(T), sys2.nx, sys1.nx)  sys2.A]
    B = [sys1.B*sys2.D ; sys2.B]
    C = [sys1.C   sys1.D*sys2.C;]
    D = [sys1.D*sys2.D;]

    return StateSpace{T,MT}(A, B, C, D, sys2.Ts)
end

*(sys::StateSpace, n::Number) = StateSpace(sys.A, sys.B, sys.C*n, sys.D*n, sys.Ts)
*(n::Number, sys::StateSpace) = *(sys, n)

## DIVISION ##
/(sys1::StateSpace, sys2::StateSpace) = sys1*inv(sys2)

function /(n::Number, sys::StateSpace)
    # Ensure s.D is invertible
    A, B, C, D = ssdata(sys)
    Dinv = try
        inv(D)
    catch
        error("D isn't invertible")
    end
    return StateSpace(A - B*Dinv*C, B*Dinv, -n*Dinv*C, n*Dinv, get_Ts(sys))
end

Base.inv(sys::StateSpace) = 1/sys
/(sys::StateSpace, n::Number) = StateSpace(sys.A, sys.B, sys.C/n, sys.D/n, sys.Ts)

Base.:^(sys::StateSpace, p::Integer) = Base.power_by_squaring(sys, p)


#####################################################################
##                       Indexing Functions                        ##
#####################################################################
Base.ndims(::StateSpace) = 2 # NOTE: Also for SISO systems?
Base.size(sys::StateSpace) = (noutputs(sys), ninputs(sys)) # NOTE: or just size(get_D(sys))
Base.size(sys::StateSpace, d) = d <= 2 ? size(sys)[d] : 1
Base.eltype(::Type{S}) where {S<:StateSpace} = S

function Base.getindex(sys::StateSpace, inds...)
    if size(inds, 1) != 2
        error("Must specify 2 indices to index statespace model")
    end
    rows, cols = index2range(inds...) # FIXME: ControlSystems.index2range(inds...)
    return StateSpace(copy(sys.A), sys.B[:, cols], sys.C[rows, :], sys.D[rows, cols], sys.Ts)
end

#####################################################################
##                        Display Functions                        ##
#####################################################################

function _string_mat_with_headers(X::Matrix)
    #mat = [[""] reshape(cols,1,length(cols));
    #       rows X]
    p = (io, m) -> Base.print_matrix(io, m)
    return replace(sprint(p, X), "\"" => " ")
end

Base.print(io::IO, sys::StateSpace) = show(io, sys)

function Base.show(io::IO, sys::StateSpace)
    # Compose the name vectors
    #inputs = format_names(s.inputnames, "u", "?")
    #outputs = format_names(s.outputnames, "y", "?")
    #println(io, "StateSpace:")
    println(io, typeof(sys))
    if nstates(sys) > 0
        #states = format_names(s.statenames, "x", "?")
        println(io, "A = \n", _string_mat_with_headers(sys.A))
        println(io, "B = \n", _string_mat_with_headers(sys.B))
        println(io, "C = \n", _string_mat_with_headers(sys.C))
    end
    println(io, "D = \n", _string_mat_with_headers(sys.D), "\n")
    # Print sample time
    if sys.Ts > 0
        println(io, "Sample Time: ", sys.Ts, " (seconds)")
    elseif sys.Ts == -1
        println(io, "Sample Time: unspecified")
    end
    # Print model type
    if nstates(sys) == 0
        print(io, "Static gain") # NOTE: Not quite...still has a time type
    elseif iscontinuous(sys)
        print(io, "Continuous-time state-space model")
    else
        print(io, "Discrete-time state-space model")
    end
end




"""
`minsys = minreal(s::StateSpace, tol=sqrt(eps()))` is implemented via `baltrunc` and returns a system on diagonal form.
"""
function minreal(s::StateSpace, tol=sqrt(eps()))
    s = baltrunc(s, atol=tol, rtol = 0)[1]
    try
        return diagonalize(s)
    catch
        error("Minreal only implemented for diagonalizable systems.")
    end
end



"""
`dsys = diagonalize(s::StateSpace, digits=12)` Diagonalizes the system such that the A-matrix is diagonal.
"""
function diagonalize(s::StateSpace, digits = 12)
    r = x -> round(x, digits=digits)
    S,V = eigen(s.A)
    try
        A = diagm(0 => S)
        B = V\s.B
        C = s.C*V
        D = s.D
        return ss(A,B,C,D)
    catch e
        error("System not diagonalizable", e)
    end
end
