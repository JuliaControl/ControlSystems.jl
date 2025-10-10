struct BodemagWorkspace{T}
    R::Array{Complex{T}, 3}
    mag::Array{T, 3}
end
function BodemagWorkspace{T}(ny::Int, nu::Int, N::Int) where T <: Real
    R = Array{Complex{T},3}(undef, ny, nu, N)
    mag = Array{T,3}(undef, ny, nu, N)
    BodemagWorkspace{T}(R, mag)
end

"""
    BodemagWorkspace(sys::LTISystem, N::Int)
    BodemagWorkspace(sys::LTISystem, ω::AbstractVector)
    BodemagWorkspace(R::Array{Complex{T}, 3}, mag::Array{T, 3})
    BodemagWorkspace{T}(ny, nu, N)

Generate a workspace object for use with the in-place function [`bodemag!`](@ref).
`N` is the number of frequency points, alternatively, the input `ω` can be provided instead of `N`.
Note: for threaded applications, create one workspace object per thread. 

# Arguments:
- `mag`: The output array ∈ 𝐑(ny, nu, nω)
- `R`: Frequency-response array ∈ 𝐂(ny, nu, nω)
"""
function BodemagWorkspace(sys::LTISystem, N::Int)
    T = float(numeric_type(sys))
    R = Array{Complex{T},3}(undef, sys.ny, sys.nu, N)
    mag = Array{T,3}(undef, sys.ny, sys.nu, N)
    BodemagWorkspace(R, mag)
end

BodemagWorkspace(sys::LTISystem, ω::AbstractVector) = BodemagWorkspace(sys, length(ω))

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

Evaluate the frequency response of a linear system.

For continuous systems, computes `G(jω) = C(jωI - A)^{-1}B + D`
For discrete systems, computes `G(e^{jωT}) = C(e^{jωT}I - A)^{-1}B + D`

# Arguments
- `sys::LTISystem`: The system to analyze
- `w::AbstractVector{<:Real}`: Frequency vector (rad/s)

# Returns
- `sys_fr`: Complex frequency response array of size `(ny, nu, length(w))`

See also [`freqresp!`](@ref), [`evalfr`](@ref), [`bode`](@ref), [`nyquist`](@ref).

# Examples
```julia
using ControlSystemsBase
sys = ss([-1 0; 0 -2], [1; 1], [1 1], 0)
w = exp10.(LinRange(-2, 2, 200))
resp = freqresp(sys, w)
```
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
    R
end

function freqresp!(R::Array{T,3}, sys::TransferFunction, w_vec::AbstractVector{<:Real}) where T
    te = sys.timeevol
    ny,nu = noutputs(sys), ninputs(sys)
    @boundscheck size(R) == (ny,nu,length(w_vec))
    @inbounds for wi = eachindex(w_vec), ui = 1:nu, yi = 1:ny
        R[yi,ui,wi] = evalfr(sys.matrix[yi,ui], _freq(w_vec[wi], te))
    end
    R
end

@autovec () function freqresp(G::AbstractMatrix, w_vec::AbstractVector{<:Real})
    repeat(G, 1, 1, length(w_vec))
end

@autovec () function freqresp(G::Number, w_vec::AbstractVector{<:Real})
    fill(G, 1, 1, length(w_vec))
end

_freq(w, ::Continuous) = complex(0, w)
_freq(w, te::Discrete) = cis(w*te.Ts)

@autovec () function freqresp!(R::Array{T,3}, sys::AbstractStateSpace, w_vec::AbstractVector{W}) where {T, W <: Real}
    ny, nu = size(sys)
    @boundscheck size(R) == (ny,nu,length(w_vec))
    if sys.nx == 0 # Only D-matrix
        @inbounds for i in eachindex(w_vec)
            R[:,:,i] .= sys.D
        end
        return R
    end
    local F, Q
    try
        F = hessenberg(sys.A)
        Q = Matrix(F.Q)
    catch e
        # For matrix types that do not have a hessenberg implementation, we call the standard version of freqresp.
        (e isa @static VERSION < v"1.12" ? Union{MethodError, ErrorException} : Union{MethodError, ErrorException, FieldError}) && return freqresp_nohess!(R, sys, w_vec)
        # ErrorException appears if we try to access Q on a type which does not have Q as a field or property, notably HessenbergFactorization from GenericLinearAlgebra, on julia v1.12, this is instead a FieldError
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
    R
end

#=
Custom implementation of hessenberg ldiv to allow reuse of u and cs
With this method, the following benchmark goes from 
(100017 allocations: 11.48 MiB) # before
to 
(17 allocations: 35.45 KiB)     # after

w = exp10.(LinRange(-2, 2, 50000))
G = ssrand(2,2,3)
R = freqresp(G, w);
@btime ControlSystemsBase.freqresp!($R, $G, $w);
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
    if nx == 0 # Only D-matrix
        @inbounds for i in eachindex(w_vec)
            R[:,:,i] .= sys.D
        end
        return R
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
        # @views copyto!(Ac[Adiag],A[Adiag]) # reset storage to A
        for j in Adiag # Loop slightly faster
            Ac[j] = A[j] # reset storage to A
        end
        @views Ac[Adiag] .-= w # Ac = A - w*I
        Bc = Ac \ B # Bc = (A - w*I)\B # avoid inplace to handle sparse matrices etc.
        mul!(Ri, C, Bc, -1, 1) # use of 5-arg mul to subtract from D already in Ri. - rather than + since (A - w*I) instead of (w*I - A)
    end
    R
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
    [v[i,j]  for i in 1:nu, j in 1:ny, v in vals]
end

"""
    mag, phase, w = bode(sys[, w]; unwrap=true)

Compute the magnitude and phase parts of the frequency response of system `sys`
at frequencies `w`. The frequency response is evaluated as `G(jω)` for continuous
systems and `G(e^{jωT})` for discrete systems.

# Arguments
- `sys::LTISystem`: The system to analyze
- `w::AbstractVector`: Frequency vector (rad/s). If omitted, a default frequency range is used.
- `unwrap::Bool`: If true (default), apply phase unwrapping to avoid discontinuities

# Returns
- `mag`: Magnitude of frequency response, size `(ny, nu, length(w))`
- `phase`: Phase in degrees, size `(ny, nu, length(w))`
- `w`: Frequency vector used

Note: Phase unwrapping is computationally expensive and can be disabled for better performance. For higher performance, see the function [`bodemag!`](@ref) that computes the magnitude only.

See also [`bodeplot`](@ref), [`bodemag!`](@ref), [`freqresp`](@ref).

# Examples
```julia
using ControlSystemsBase
sys = tf(1, [1, 1])
mag, phase, w = bode(sys)
```
""" 
@autovec (1, 2) function bode(sys::LTISystem, w::AbstractVector; unwrap=true)
    resp = freqresp(sys, w)
    angles = angle.(resp)
    unwrap && unwrap!(angles,3)
    @. angles = rad2deg(angles)
    return abs.(resp), angles, w
end
@autovec (1, 2) bode(sys::LTISystem) = bode(sys, _default_freq_vector(sys, Val{:bode}()))

# Performance difference between bode and bodemag for tf. Note how expensive the phase unwrapping is.
# using ControlSystemsBase
# G = tf(ssrand(2,2,5))
# w = exp10.(LinRange(-2, 2, 20000))
# @btime bode($G, $w);
# # 55.120 ms (517957 allocations: 24.42 MiB)
# @btime bode($G, $w, unwrap=false);
# # 3.624 ms (7 allocations: 2.44 MiB)
# ws = ControlSystemsBase.BodemagWorkspace(G, w)
# @btime bodemag!($ws, $G, $w);
# # 2.991 ms (1 allocation: 64 bytes)

"""
    mag = bodemag!(ws::BodemagWorkspace, sys::LTISystem, w::AbstractVector)

Compute the Bode magnitude operating in-place on an instance of [`BodemagWorkspace`](@ref).

# Arguments
- `ws::BodemagWorkspace`: Pre-allocated workspace created with [`BodemagWorkspace`](@ref)
- `sys::LTISystem`: The system to analyze
- `w::AbstractVector`: Frequency vector (rad/s)

# Returns
- `mag`: Magnitude of frequency response, size `(ny, nu, length(w))`. 
  Note: The returned array is aliased with `ws.mag`.

# Performance
This function provides significant performance benefits for repeated magnitude calculations:
- Avoids memory allocation by reusing workspace arrays
- Optimized for systems with many frequency points

For thread-safe applications, create one workspace per thread.

See also [`BodemagWorkspace`](@ref), [`bode`](@ref), [`freqresp!`](@ref).

# Examples
```julia
using ControlSystemsBase
sys = tf(1, [1, 1])
w = exp10.(LinRange(-2, 2, 200))
ws = BodemagWorkspace(sys, w)
mag = bodemag!(ws, sys, w)
```
"""
function bodemag!(ws::BodemagWorkspace, sys::LTISystem, w::AbstractVector)
    freqresp!(ws.R, sys, w)
    @. ws.mag = abs(ws.R)
    ws.mag
end

function bodemag_nohess!(ws::BodemagWorkspace, sys::LTISystem, w::AbstractVector)
    freqresp_nohess!(ws.R, sys, w)
    @. ws.mag = abs(ws.R)
    ws.mag
end

"""
    re, img, w = nyquist(sys[, w])

Compute the real and imaginary parts of the frequency response of system `sys`
at frequencies `w`. The frequency response is evaluated as `G(jω)` for continuous
systems and `G(e^{jωT})` for discrete systems.

# Arguments
- `sys::LTISystem`: The system to analyze
- `w::AbstractVector`: Frequency vector (rad/s). If omitted, a default frequency range is used.

# Returns
- `re`: Real part of frequency response, size `(ny, nu, length(w))`
- `img`: Imaginary part of frequency response, size `(ny, nu, length(w))`
- `w`: Frequency vector used

See also [`nyquistplot`](@ref), [`freqresp`](@ref), [`bode`](@ref).

# Examples
```julia
using ControlSystems
sys = tf(1, [1, 1])
w = logspace(-2, 2, 100)
re, img, w = nyquist(sys, w)
```
""" 
@autovec (1, 2) function nyquist(sys::LTISystem, w::AbstractVector)
    resp = freqresp(sys, w)
    return real(resp), imag(resp), w
end
@autovec (1, 2) nyquist(sys::LTISystem) = nyquist(sys, _default_freq_vector(sys, Val{:nyquist}()))

"""
    sv, w = sigma(sys[, w])

Compute the singular values of the frequency response of system `sys` at
frequencies `w`. The frequency response is evaluated as `G(jω)` for continuous
systems and `G(e^{jωT})` for discrete systems.

# Arguments
- `sys::LTISystem`: The system to analyze
- `w::AbstractVector`: Frequency vector (rad/s). If omitted, a default frequency range is used.

# Returns
- `sv`: Singular values of frequency response, size `(min(ny, nu), length(w))`
- `w`: Frequency vector used

The singular values provide information about the system's gain characteristics
across different frequencies and input/output directions.

See also [`sigmaplot`](@ref), [`freqresp`](@ref), [`bode`](@ref).

# Examples
```julia
using ControlSystems
sys = ss([-1 0; 0 -2], [1 0; 0 1], [1 1; 0 1], 0)
sv, w = sigma(sys)
```
""" 
@autovec (1) function sigma(sys::LTISystem, w::AbstractVector)
    resp = freqresp(sys, w)
    ny, nu = size(sys)
    if ny == 1 || nu == 1 # Shortcut available
        sv = Matrix{real(eltype(resp))}(undef, 1, length(w))
        for i = eachindex(w)
            @views sv[1, i] = norm(resp[:,:,i])
        end
    else
        sv = dropdims(mapslices(svdvals, resp, dims=(1,2)),dims=2)::Matrix{real(eltype(resp))}
    end
    return sv, w
end
@autovec (1) sigma(sys::LTISystem) = sigma(sys, _default_freq_vector(sys, Val{:sigma}()))

function _default_freq_vector(systems::Vector{<:LTISystem}, plot; adaptive=false)
    if adaptive
        min_pt_per_dec = 100
        min_pt_total = 1000
    else
        min_pt_per_dec = 60
        min_pt_total = 200
    end
    bounds = map(sys -> _bounds_and_features(sys, plot)[1], systems)
    w1 = minimum(minimum, bounds)
    w2 = maximum(maximum, bounds)

    nw = round(Int, max(min_pt_total, min_pt_per_dec*(w2 - w1)))
    w = exp10.(range(w1, stop=w2, length=nw))
    if length(systems) == 1 && isdiscrete(systems[1])
        w[end] = π/systems[1].Ts # To account for numerical rounding problems from exp(log())
    end
    w
end
_default_freq_vector(sys::LTISystem, plot; kwargs...) = _default_freq_vector(
        [sys], plot; kwargs...)

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
    if isdiscrete(sys)
        w2 = log10(π/sys.Ts) # Draw up to Nyquist frequency for discrete systems
    end
    return [w1, w2], zp
end
