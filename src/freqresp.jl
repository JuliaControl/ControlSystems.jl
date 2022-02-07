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

`w -> C*((iw*im -A)^-1)*B + D`

of system `sys` over the frequency vector `w`.
"""
@autovec () function freqresp(sys::LTISystem, w_vec::AbstractVector{<:Real})
    te = sys.timeevol
    ny,nu = noutputs(sys), ninputs(sys)
    [evalfr(sys[i,j], _freq(w, te))[] for w in w_vec, i in 1:ny, j in 1:nu]
end

@autovec () function freqresp(G::AbstractMatrix, w_vec::AbstractVector{<:Real})
    repeat(G, 1, 1, length(w_vec))
end

@autovec () function freqresp(G::Number, w_vec::AbstractVector{<:Real})
    fill(G, 1, 1, length(w_vec))
end

_freq(w, ::Continuous) = complex(0, w)
_freq(w, te::Discrete) = cis(w*te.Ts)

@autovec () function freqresp(sys::AbstractStateSpace, w_vec::AbstractVector{W}) where W <: Real
    ny, nu = size(sys)
    T = promote_type(Complex{real(eltype(sys.A))}, Complex{W})
    if sys.nx == 0 # Only D-matrix
        return PermutedDimsArray(repeat(T.(sys.D), 1, 1, length(w_vec)), (3,1,2))
    end
    local F
    try
        F = hessenberg(sys.A)
    catch e
        # For matrix types that do not have a hessenberg implementation, we call the standard version of freqresp.
        e isa MethodError && return freqresp_nohess(sys, w_vec)
        rethrow()
    end
    Q = Matrix(F.Q)
    A = F.H
    C = sys.C*Q
    B = Q\sys.B 
    D = sys.D

    te = sys.timeevol
    R = Array{T, 3}(undef, ny, nu, length(w_vec))
    Bc = similar(B, T) # for storage
    for i in eachindex(w_vec)
        Ri = @views R[:,:,i]
        copyto!(Ri,D) # start with the D-matrix
        isinf(w_vec[i]) && continue
        copyto!(Bc,B) # initialize storage to B
        w = -_freq(w_vec[i], te)
        ldiv!(A, Bc, shift = w) # B += (A - w*I)\B # solve (A-wI)X = B, storing result in B
        mul!(Ri, C, Bc, -1, 1) # use of 5-arg mul to subtract from D already in Ri. - rather than + since (A - w*I) instead of (w*I - A)
    end
    PermutedDimsArray(R, (3,1,2)) # PermutedDimsArray doesn't allocate to perform the permutation
end

"""
    freqresp_nohess(sys::AbstractStateSpace, w_vec::AbstractVector{<:Real})

Compute the frequency response of `sys` without forming a Hessenberg factorization.
This function is called automatically if the Hessenberg factorization fails.
"""
freqresp_nohess
@autovec () function freqresp_nohess(sys::AbstractStateSpace, w_vec::AbstractVector{W}) where W <: Real
    ny, nu = size(sys)
    nx = sys.nx
    T = promote_type(Complex{real(eltype(sys.A))}, Complex{W})
    if nx == 0 # Only D-matrix
        return PermutedDimsArray(repeat(T.(sys.D), 1, 1, length(w_vec)), (3,1,2))
    end
    A,B,C,D = ssdata(sys)
    te = sys.timeevol
    R = Array{T, 3}(undef, ny, nu, length(w_vec))
    Ac = (A+one(T)*I) # for storage
    Adiag = diagind(A)
    for i in eachindex(w_vec)
        Ri = @views R[:,:,i]
        copyto!(Ri,D) # start with the D-matrix
        isinf(w_vec[i]) && continue
        w = _freq(w_vec[i], te)
        @views copyto!(Ac[Adiag],A[Adiag]) # reset storage to A
        @views Ac[Adiag] .-= w # Ac = A - w*I
        Bc = Ac \ B # Bc = (A - w*I)\B # avoid inplace to handle sparse matrices etc.
        mul!(Ri, C, Bc, -1, 1) # use of 5-arg mul to subtract from D already in Ri. - rather than + since (A - w*I) instead of (w*I - A)
    end
    PermutedDimsArray(R, (3,1,2)) # PermutedDimsArray doesn't allocate to perform the permutation
end


function _evalfr_return_type(sys::AbstractStateSpace, s::Number)
    T0 = numeric_type(sys)
    temp_product = one(T0)*one(typeof(s))
    typeof(temp_product/temp_product)
end

"""
`evalfr(sys, x)` Evaluate the transfer function of the LTI system sys
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

"""`mag, phase, w = bode(sys[, w])`

Compute the magnitude and phase parts of the frequency response of system `sys`
at frequencies `w`

`mag` and `phase` has size `(length(w), ny, nu)`""" 
@autovec (1, 2) function bode(sys::LTISystem, w::AbstractVector)
    resp = freqresp(sys, w)
    return abs.(resp), rad2deg.(unwrap!(angle.(resp),1)), w
end
@autovec (1, 2) bode(sys::LTISystem) = bode(sys, _default_freq_vector(sys, Val{:bode}()))

"""`re, im, w = nyquist(sys[, w])`

Compute the real and imaginary parts of the frequency response of system `sys`
at frequencies `w`

`re` and `im` has size `(length(w), ny, nu)`""" 
@autovec (1, 2) function nyquist(sys::LTISystem, w::AbstractVector)
    resp = freqresp(sys, w)
    return real(resp), imag(resp), w
end
@autovec (1, 2) nyquist(sys::LTISystem) = nyquist(sys, _default_freq_vector(sys, Val{:nyquist}()))

"""`sv, w = sigma(sys[, w])`

Compute the singular values `sv` of the frequency response of system `sys` at
frequencies `w`

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
