function freqresp(sys::LTISystem, w::Real)
    # Create imaginary freq vector s
    if iscontinuous(sys)
        s = im*w
    else
        s = cis(w*sys.Ts)
    end
    evalfr(sys, s)
end

freqresp(G::Union{UniformScaling, AbstractMatrix, Number}, w::Real) = G

"""
    sys_fr = freqresp(sys, w)

Evaluate the frequency response of a linear system

`w -> C*((iw*im*I - A)^-1)*B + D`

of system `sys` over the frequency vector `w`.
"""
@autovec () function freqresp(sys::LTISystem, w_vec::AbstractVector{W}) where W <: Real
    te = timeevol(sys)
    ny,nu = noutputs(sys), ninputs(sys)
    T = promote_type(Complex{real(numeric_type(sys))}, Complex{W})
    R = Array{T, 3}(undef, ny, nu, length(w_vec))
    freqresp!(R, sys, w_vec)
end

"""
    freqresp!(R::Array{T, 3}, sys::LTISystem, w_vec::AbstractVector{<:Real})

In-place version of [`freqresp`](@ref) that takes a pre-allocated array `R` of size (ny, nu, nw)`
"""
function freqresp!(R::Array{T,3}, sys::LTISystem, w_vec::AbstractVector{<:Real}) where T
    te = sys.timeevol
    ny,nu = noutputs(sys), ninputs(sys)
    @boundscheck size(R) == (ny,nu,length(w_vec))
    @inbounds for wi = eachindex(w_vec), ui = 1:nu, yi = 1:ny
        R[yi,ui,wi] = evalfr(sys[yi,ui], _freq(w_vec[wi], te))[]
    end
    PermutedDimsArray{T,3,(3,1,2),(2,3,1),Array{T,3}}(R)
end

@autovec () function freqresp(G::AbstractMatrix, w_vec::AbstractVector{<:Real})
    repeat(G, 1, 1, length(w_vec))
end

@autovec () function freqresp(G::Number, w_vec::AbstractVector{<:Real})
    fill(G, 1, 1, length(w_vec))
end

_freq(w, ::Continuous) = complex(0, w)
_freq(w, te::Discrete) = cis(w*te.Ts)

function freqresp(sys::AbstractStateSpace, w_vec::AbstractVector{W}) where W <: Real
    ny, nu = size(sys)
    T = promote_type(Complex{real(eltype(sys.A))}, Complex{W})
    R = Array{T, 3}(undef, ny, nu, length(w_vec))
    freqresp!(R, sys, w_vec)
end
@autovec () function freqresp!(R::Array{T,3}, sys::AbstractStateSpace, w_vec::AbstractVector{W}) where {T, W <: Real}
    ny, nu = size(sys)
    @boundscheck size(R) == (ny,nu,length(w_vec))
    PDT = PermutedDimsArray{T,3,(3,1,2),(2,3,1),Array{T,3}}
    if sys.nx == 0 # Only D-matrix
        @inbounds for i in eachindex(w_vec)
            R[:,:,i] .= sys.D
        end
        return PDT(R)::PDT
    end
    local F, Q
    try
        F = hessenberg(sys.A)
        Q = Matrix(F.Q)
    catch e
        # For matrix types that do not have a hessenberg implementation, we call the standard version of freqresp.
        e isa Union{MethodError, ErrorException} && return freqresp_nohess!(R, sys, w_vec)
        # ErrorException appears if we try to access Q on a type which does not have Q as a field or property, notably HessenbergFactorization from GenericLinearAlgebra
        rethrow()
    end
    A = F.H
    C = complex.(sys.C*Q) # We make C complex in order to not incur allocations in mul! below
    B = Q\sys.B 
    D = sys.D

    te = sys.timeevol
    Bc = similar(B, T) # for storage
    w = -_freq(w_vec[1], te)
    u = Vector{typeof(zero(eltype(A.data))+w)}(undef, sys.nx)
    cs = Vector{Tuple{real(eltype(u)),eltype(u)}}(undef, length(u)) # store Givens rotations
    @inbounds for i in eachindex(w_vec)
        Ri = @view(R[:, :, i])
        copyto!(Ri,D) # start with the D-matrix
        isinf(w_vec[i]) && continue
        copyto!(Bc,B) # initialize storage to B
        w = -_freq(w_vec[i], te)
        ldiv2!(u, cs, A, Bc, shift = w) # B += (A - w*I)\B # solve (A-wI)X = B, storing result in B
        mul!(Ri, C, Bc, -1, 1) # use of 5-arg mul to subtract from D already in Ri. - rather than + since (A - w*I) instead of (w*I - A)
    end
    PDT(R)::PDT # PermutedDimsArray doesn't allocate to perform the permutation
end

#=
Custom implementation of hessenberg ldiv to allow reuse of u and cs
With this method, the following benchmark goes from 
(100017 allocations: 11.48 MiB) # before
to 
(17 allocations: 35.45 KiB)     # after

w = exp10.(LinRange(-2, 2, 50000))
G = ssrand(2,2,3)
R = freqresp(G, w).parent;
@btime ControlSystems.freqresp!($R, $G, $w);
=# 
function ldiv2!(u, cs, F::UpperHessenberg, B::AbstractVecOrMat; shift::Number=false)
    LinearAlgebra.checksquare(F)
    m = size(F,1)
    m != size(B,1) && throw(DimensionMismatch("wrong right-hand-side # rows != $m"))
    LinearAlgebra.require_one_based_indexing(B)
    n = size(B,2)
    H = F.data
    μ = shift
    copyto!(u, 1, H, m*(m-1)+1, m) # u .= H[:,m]
    u[m] += μ
    X = B # not a copy, just rename to match paper
    @inbounds for k = m:-1:2
        c, s, ρ = LinearAlgebra.givensAlgorithm(u[k], H[k,k-1])
        cs[k] = (c, s)
        for i = 1:n
            X[k,i] /= ρ
            t₁ = s * X[k,i]; t₂ = c * X[k,i]
            @simd for j = 1:k-2
                X[j,i] -= u[j]*t₂ + H[j,k-1]*t₁
            end
            X[k-1,i] -= u[k-1]*t₂ + (H[k-1,k-1] + μ) * t₁
        end
        @simd for j = 1:k-2
            u[j] = H[j,k-1]*c - u[j]*s'
        end
        u[k-1] = (H[k-1,k-1] + μ) * c - u[k-1]*s'
    end
    for i = 1:n
        τ₁ = X[1,i] / u[1]
        @inbounds for j = 2:m
            τ₂ = X[j,i]
            c, s = cs[j]
            X[j-1,i] = c*τ₁ + s*τ₂
            τ₁ = c*τ₂ - s'τ₁
        end
        X[m,i] = τ₁
    end
    return X
end

function freqresp_nohess(sys::AbstractStateSpace, w_vec::AbstractVector{W}) where W <: Real
    ny, nu = size(sys)
    T = promote_type(Complex{real(eltype(sys.A))}, Complex{W})
    R = Array{T, 3}(undef, ny, nu, length(w_vec))
    freqresp_nohess!(R, sys, w_vec)
end

"""
    freqresp_nohess(sys::AbstractStateSpace, w_vec::AbstractVector{<:Real})

Compute the frequency response of `sys` without forming a Hessenberg factorization.
This function is called automatically if the Hessenberg factorization fails.
"""
freqresp_nohess
@autovec () function freqresp_nohess!(R::Array{T,3}, sys::AbstractStateSpace, w_vec::AbstractVector{W}) where {T, W <: Real}
    ny, nu = size(sys)
    @boundscheck size(R) == (ny,nu,length(w_vec))
    nx = sys.nx
    PDT = PermutedDimsArray{T,3,(3,1,2),(2,3,1),Array{T,3}}
    if nx == 0 # Only D-matrix
        @inbounds for i in eachindex(w_vec)
            R[:,:,i] .= sys.D
        end
        return PDT(R)::PDT
    end
    A,B,C0,D = ssdata(sys)
    C = complex.(C0) # We make C complex in order to not incur allocations in mul! below
    te = sys.timeevol
    Ac = (A+one(T)*I) # for storage
    Adiag = diagind(A)
    @inbounds for i in eachindex(w_vec)
        Ri = @views R[:,:,i]
        copyto!(Ri,D) # start with the D-matrix
        isinf(w_vec[i]) && continue
        w = _freq(w_vec[i], te)
        @views copyto!(Ac[Adiag],A[Adiag]) # reset storage to A
        @views Ac[Adiag] .-= w # Ac = A - w*I
        Bc = Ac \ B # Bc = (A - w*I)\B # avoid inplace to handle sparse matrices etc.
        mul!(Ri, C, Bc, -1, 1) # use of 5-arg mul to subtract from D already in Ri. - rather than + since (A - w*I) instead of (w*I - A)
    end
    PDT(R)::PDT # PermutedDimsArray doesn't allocate to perform the permutation
end


function _evalfr_return_type(sys::AbstractStateSpace, s::Number)
    T0 = numeric_type(sys)
    temp_product = one(T0)*one(typeof(s))
    typeof(temp_product/temp_product)
end

"""
    evalfr(sys, x)
    
Evaluate the transfer function of the LTI system sys
at the complex number s=x (continuous-time) or z=x (discrete-time).

For many values of `x`, use `freqresp` instead.
"""
function evalfr(sys::AbstractStateSpace, s::Number)
    T = _evalfr_return_type(sys, s)
    try
        R = s*I - sys.A
        sys.D + sys.C*((R\sys.B))
    catch e
        @warn "Got exception $e, returning Inf" maxlog=1
        fill(convert(T, Inf), size(sys))
    end
end

function evalfr(G::TransferFunction{<:TimeEvolution,<:SisoTf}, s::Number)
    map(m -> evalfr(m,s), G.matrix)
end

evalfr(G::Union{UniformScaling, AbstractMatrix, Number}, s) = G

"""
`F(s)`, `F(omega, true)`, `F(z, false)`

Notation for frequency response evaluation.
- F(s) evaluates the continuous-time transfer function F at s.
- F(omega,true) evaluates the discrete-time transfer function F at exp(im*Ts*omega)
- F(z,false) evaluates the discrete-time transfer function F at z
"""
function (sys::TransferFunction)(s)
    evalfr(sys,s)
end

function (sys::TransferFunction)(z_or_omega::Number, map_to_unit_circle::Bool)
    isdiscrete(sys) || throw(ArgumentError("It only makes no sense to call this function with discrete systems"))
    if map_to_unit_circle
        isreal(z_or_omega) ? evalfr(sys,exp(im*z_or_omega.*sys.Ts)) : error("To map to the unit circle, omega should be real")
    else
        evalfr(sys,z_or_omega)
    end
end

function (sys::TransferFunction)(z_or_omegas::AbstractVector, map_to_unit_circle::Bool)
    isdiscrete(sys) || throw(ArgumentError("It only makes no sense to call this function with discrete systems"))
    vals = sys.(z_or_omegas, map_to_unit_circle)# evalfr.(sys,exp.(evalpoints))
    # Reshape from vector of evalfr matrizes, to (in,out,freq) Array
    nu,ny = size(vals[1])
    [v[i,j]  for v in vals, i in 1:nu, j in 1:ny]
end

"""
    mag, phase, w = bode(sys[, w])

Compute the magnitude and phase parts of the frequency response of system `sys`
at frequencies `w`. See also [`bodeplot`](@ref)

`mag` and `phase` has size `(length(w), ny, nu)`""" 
@autovec (1, 2) function bode(sys::LTISystem, w::AbstractVector)
    resp = freqresp(sys, w)
    return abs.(resp), rad2deg.(unwrap!(angle.(resp),1)), w
end
@autovec (1, 2) bode(sys::LTISystem) = bode(sys, _default_freq_vector(sys, Val{:bode}()))

"""
    re, im, w = nyquist(sys[, w])

Compute the real and imaginary parts of the frequency response of system `sys`
at frequencies `w`. See also [`nyquistplot`](@ref)

`re` and `im` has size `(length(w), ny, nu)`""" 
@autovec (1, 2) function nyquist(sys::LTISystem, w::AbstractVector)
    resp = freqresp(sys, w)
    return real(resp), imag(resp), w
end
@autovec (1, 2) nyquist(sys::LTISystem) = nyquist(sys, _default_freq_vector(sys, Val{:nyquist}()))

"""
    sv, w = sigma(sys[, w])

Compute the singular values `sv` of the frequency response of system `sys` at
frequencies `w`. See also [`sigmaplot`](@ref)

`sv` has size `(length(w), max(ny, nu))`""" 
@autovec (1) function sigma(sys::LTISystem, w::AbstractVector)
    resp = freqresp(sys, w)
    sv = dropdims(mapslices(svdvals, resp, dims=(2,3)),dims=3)
    return sv, w
end
@autovec (1) sigma(sys::LTISystem) = sigma(sys, _default_freq_vector(sys, Val{:sigma}()))

function _default_freq_vector(systems::Vector{<:LTISystem}, plot)
    min_pt_per_dec = 60
    min_pt_total = 200
    bounds = map(sys -> _bounds_and_features(sys, plot)[1], systems)
    w1 = minimum(minimum.(bounds))
    w2 = maximum(maximum.(bounds))

    nw = round(Int, max(min_pt_total, min_pt_per_dec*(w2 - w1)))
    return exp10.(range(w1, stop=w2, length=nw))
end
_default_freq_vector(sys::LTISystem, plot) = _default_freq_vector(
        [sys], plot)


function _bounds_and_features(sys::LTISystem, plot::Val)
    # Get zeros and poles for each channel
    if !isa(plot, Val{:sigma})
        zs, ps = zpkdata(sys)
        # Compose vector of all zs, ps, positive conjugates only.
        zpType = promote_type(eltype(eltype(zs)), eltype(eltype(ps)))
        zp = vcat(zpType[], zs..., ps...) # Emty vector to avoid type unstable vcat()
        zp = zp[imag(zp) .>= 0.0]
    else
        # For sigma plots, use the MIMO poles and zeros
        zp = [tzeros(sys); poles(sys)]
    end
    # Get the frequencies of the features, ignoring low frequency dynamics
    fzp = log10.(abs.(zp))
    fzp = fzp[fzp .> -4]
    fzp = sort!(fzp)
    # Determine the bounds on the frequency vector
    if !isempty(fzp)
        w1 = floor(fzp[1] - 1.2)
        w2 = ceil(fzp[end] + 1.2)
        # Expand the range for nyquist plots
        if plot isa Val{:nyquist}
            w1 -= 0.0
            w2 += 1.0
        end
    else
        w1 = 0.0
        w2 = 2.0
    end
    if isdiscrete(sys) # Do not draw above Nyquist freq for disc. systems
        w2 = min(w2, log10(π/sys.Ts))
    end
    return [w1, w2], zp
end
