"""sys_fr = freqresp(sys, w)

Evaluate the frequency response of a linear system

`w -> C*((iw*im -A)^-1)*B + D`

of system `sys` over the frequency vector `w`."""
@autovec () function freqresp(sys::LTISystem, w_vec::AbstractVector{<:Real})
    # Create imaginary freq vector s
    if iscontinuous(sys)
        s_vec = im*w_vec
    else
        s_vec = exp.(w_vec*(im*sys.Ts))
    end
    if isa(sys, StateSpace)
        sys = _preprocess_for_freqresp(sys)
    end
    ny,nu = noutputs(sys), ninputs(sys)
    [evalfr(sys[i,j], s)[] for s in s_vec, i in 1:ny, j in 1:nu]
end

# Implements algorithm found in:
# Laub, A.J., "Efficient Multivariable Frequency Response Computations",
# IEEE Transactions on Automatic Control, AC-26 (1981), pp. 407-408.
function _preprocess_for_freqresp(sys::StateSpace)
    if isempty(sys.A) # hessfact does not work for empty matrices
        return sys
    end
    Tsys = numeric_type(sys)
    TT = promote_type(typeof(zero(Tsys)/norm(one(Tsys))), Float32)

    A, B, C, D = sys.A, sys.B, sys.C, sys.D
    F = hessenberg(A)
    T = F.Q
    P = C*T
    Q = T\B # TODO Type stability? # T is unitary, so mutliplication with T' should do the trick
    # FIXME; No performance improvement from Hessienberg structure, also weired renaming of matrices
    StateSpace(F.H, Q, P, D, sys.timeevol)
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
        @warn "Got exception $e, returning Inf" max_log=1
        fill(convert(T, Inf), size(sys))
    end
end

function evalfr(G::TransferFunction{<:TimeEvolution,<:SisoTf}, s::Number)
    map(m -> evalfr(m,s), G.matrix)
end

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
    @assert isdiscrete(sys) "It only makes no sense to call this function with discrete systems"
    if map_to_unit_circle
        isreal(z_or_omega) ? evalfr(sys,exp(im*z_or_omega.*sys.Ts)) : error("To map to the unit circle, omega should be real")
    else
        evalfr(sys,z_or_omega)
    end
end

function (sys::TransferFunction)(z_or_omegas::AbstractVector, map_to_unit_circle::Bool)
    @assert isdiscrete(sys) "It only makes no sense to call this function with discrete systems"
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
         zp = [tzero(sys); pole(sys)]
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
            w1 -= 1.0
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
