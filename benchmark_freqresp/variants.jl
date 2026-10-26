# Candidate freqresp! implementations for benchmarking.
#
# Internals `_freq` and `ldiv2!` are unexported in ControlSystemsBase, so the modified
# variants copy them here (verbatim for the serial baseline, real/imag-split for @turbo)
# so they can be mutated freely without touching the package source.
#
# The benchmark target system is a CONTINUOUS state-space system, so the frequency map is
# `_freqc(w) = complex(0, w)` throughout. (Discrete would use `cis(w*Ts)`.)

using LinearAlgebra
using LinearAlgebra: UpperHessenberg, givensAlgorithm
import LoopVectorization                      # provides @turbo
using LoopVectorization: @turbo, @tturbo
import Polyester                              # provides @batch
import OhMyThreads

_freqc(w) = complex(0.0, w)

# ---------------------------------------------------------------------------
# Copied verbatim from lib/ControlSystemsBase/src/freqresp.jl:168-206
# (renamed to keep it independent of the package's own definition)
# ---------------------------------------------------------------------------
function ldiv2_orig!(u, cs, F::UpperHessenberg, B::AbstractVecOrMat; shift::Number=false)
    LinearAlgebra.checksquare(F)
    m = size(F,1)
    m != size(B,1) && throw(DimensionMismatch("wrong right-hand-side # rows != $m"))
    LinearAlgebra.require_one_based_indexing(B)
    n = size(B,2)
    H = F.data
    μ = shift
    copyto!(u, 1, H, m*(m-1)+1, m) # u .= H[:,m]
    u[m] += μ
    X = B
    @inbounds for k = m:-1:2
        c, s, ρ = givensAlgorithm(u[k], H[k,k-1])
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

# ---------------------------------------------------------------------------
# Real/imag-split rewrite enabling @turbo on the two hot loops.
# @turbo does NOT support Complex, so all complex arithmetic is written out by hand
# over separate real/imag Float64 arrays. H is real already; givens `c` is real, `s`/`ρ`
# complex. The two loops over `j = 1:k-2` (X update and u update) carry no cross-iteration
# dependency, so they are valid @turbo targets.
#   u, X stored split as (ur,ui) and (Xr,Xi); cs stores (c::Real, s::Complex).
# ---------------------------------------------------------------------------
function ldiv2_turbo!(ur, ui, cs, H::AbstractMatrix{<:Real}, Xr, Xi, μ::Complex)
    m = size(H, 1)
    n = size(Xr, 2)
    μr = real(μ); μi = imag(μ)
    @inbounds for j in 1:m            # u .= H[:,m]; then u[m] += μ
        ur[j] = H[j,m]; ui[j] = 0.0
    end
    ur[m] += μr; ui[m] += μi
    @inbounds for k = m:-1:2
        c, s, ρ = givensAlgorithm(complex(ur[k], ui[k]), H[k,k-1])
        cs[k] = (c, s)
        sr = real(s); si = imag(s)
        ρr = real(ρ); ρi = imag(ρ)
        ρd = ρr*ρr + ρi*ρi
        hr = H[k-1,k-1] + μr; hi = μi   # (H[k-1,k-1] + μ)
        for i = 1:n
            xr = Xr[k,i]; xi = Xi[k,i]              # X[k,i] /= ρ
            nxr = (xr*ρr + xi*ρi)/ρd
            nxi = (xi*ρr - xr*ρi)/ρd
            Xr[k,i] = nxr; Xi[k,i] = nxi
            t1r = sr*nxr - si*nxi; t1i = sr*nxi + si*nxr   # t₁ = s*X[k,i]
            t2r = c*nxr;           t2i = c*nxi             # t₂ = c*X[k,i]
            @turbo for j = 1:k-2                           # X[j,i] -= u[j]*t₂ + H[j,k-1]*t₁
                Xr[j,i] -= ur[j]*t2r - ui[j]*t2i + H[j,k-1]*t1r
                Xi[j,i] -= ur[j]*t2i + ui[j]*t2r + H[j,k-1]*t1i
            end
            # X[k-1,i] -= u[k-1]*t₂ + (H[k-1,k-1]+μ)*t₁
            Xr[k-1,i] -= ur[k-1]*t2r - ui[k-1]*t2i + (hr*t1r - hi*t1i)
            Xi[k-1,i] -= ur[k-1]*t2i + ui[k-1]*t2r + (hr*t1i + hi*t1r)
        end
        @turbo for j = 1:k-2          # u[j] = H[j,k-1]*c - u[j]*s'  (s' = conj(s))
            urj = ur[j]; uij = ui[j]
            ur[j] = H[j,k-1]*c - (urj*sr + uij*si)
            ui[j] = -(uij*sr - urj*si)
        end
        urk = ur[k-1]; uik = ui[k-1]  # u[k-1] = (H[k-1,k-1]+μ)*c - u[k-1]*s'
        ur[k-1] = hr*c - (urk*sr + uik*si)
        ui[k-1] = hi*c - (uik*sr - urk*si)
    end
    # NOTE: this back-substitution loop is NOT @turbo-able. The inner j-loop is a
    # sequential recurrence (τ₁ carried across iterations — a scan), which @turbo cannot
    # vectorize, and the only independent axis (columns i) is length 1 for SISO. Annotating
    # it with @turbo merely trips check_args and falls back to a plain @inbounds @fastmath
    # loop (no actual vectorization).
    @inbounds for i = 1:n             # back substitution
        u1r = ur[1]; u1i = ui[1]
        d = u1r*u1r + u1i*u1i
        x1r = Xr[1,i]; x1i = Xi[1,i]
        τ1r = (x1r*u1r + x1i*u1i)/d    # τ₁ = X[1,i]/u[1]
        τ1i = (x1i*u1r - x1r*u1i)/d
        for j = 2:m
            τ2r = Xr[j,i]; τ2i = Xi[j,i]
            c, s = cs[j]; sr = real(s); si = imag(s)
            Xr[j-1,i] = c*τ1r + (sr*τ2r - si*τ2i)   # X[j-1,i] = c*τ₁ + s*τ₂
            Xi[j-1,i] = c*τ1i + (sr*τ2i + si*τ2r)
            nτ1r = c*τ2r - (sr*τ1r + si*τ1i)        # τ₁ = c*τ₂ - s'τ₁
            nτ1i = c*τ2i - (sr*τ1i - si*τ1r)
            τ1r = nτ1r; τ1i = nτ1i
        end
        Xr[m,i] = τ1r; Xi[m,i] = τ1i
    end
    return
end

# ---------------------------------------------------------------------------
# Same as ldiv2_turbo! but with @tturbo (LoopVectorization's threaded @turbo) on the
# two hot j-loops. Expected to be slower here: the loops are only k-2 <= 46 long and
# sit inside a sequentially-dependent outer k loop, so per-loop thread dispatch overhead
# dominates. Measured for completeness.
# ---------------------------------------------------------------------------
function ldiv2_tturbo!(ur, ui, cs, H::AbstractMatrix{<:Real}, Xr, Xi, μ::Complex)
    m = size(H, 1)
    n = size(Xr, 2)
    μr = real(μ); μi = imag(μ)
    @inbounds for j in 1:m
        ur[j] = H[j,m]; ui[j] = 0.0
    end
    ur[m] += μr; ui[m] += μi
    @inbounds for k = m:-1:2
        c, s, ρ = givensAlgorithm(complex(ur[k], ui[k]), H[k,k-1])
        cs[k] = (c, s)
        sr = real(s); si = imag(s)
        ρr = real(ρ); ρi = imag(ρ)
        ρd = ρr*ρr + ρi*ρi
        hr = H[k-1,k-1] + μr; hi = μi
        for i = 1:n
            xr = Xr[k,i]; xi = Xi[k,i]
            nxr = (xr*ρr + xi*ρi)/ρd
            nxi = (xi*ρr - xr*ρi)/ρd
            Xr[k,i] = nxr; Xi[k,i] = nxi
            t1r = sr*nxr - si*nxi; t1i = sr*nxi + si*nxr
            t2r = c*nxr;           t2i = c*nxi
            @tturbo for j = 1:k-2
                Xr[j,i] -= ur[j]*t2r - ui[j]*t2i + H[j,k-1]*t1r
                Xi[j,i] -= ur[j]*t2i + ui[j]*t2r + H[j,k-1]*t1i
            end
            Xr[k-1,i] -= ur[k-1]*t2r - ui[k-1]*t2i + (hr*t1r - hi*t1i)
            Xi[k-1,i] -= ur[k-1]*t2i + ui[k-1]*t2r + (hr*t1i + hi*t1r)
        end
        @tturbo for j = 1:k-2
            urj = ur[j]; uij = ui[j]
            ur[j] = H[j,k-1]*c - (urj*sr + uij*si)
            ui[j] = -(uij*sr - urj*si)
        end
        urk = ur[k-1]; uik = ui[k-1]
        ur[k-1] = hr*c - (urk*sr + uik*si)
        ui[k-1] = hi*c - (uik*sr - urk*si)
    end
    @inbounds for i = 1:n
        u1r = ur[1]; u1i = ui[1]
        d = u1r*u1r + u1i*u1i
        x1r = Xr[1,i]; x1i = Xi[1,i]
        τ1r = (x1r*u1r + x1i*u1i)/d
        τ1i = (x1i*u1r - x1r*u1i)/d
        for j = 2:m
            τ2r = Xr[j,i]; τ2i = Xi[j,i]
            c, s = cs[j]; sr = real(s); si = imag(s)
            Xr[j-1,i] = c*τ1r + (sr*τ2r - si*τ2i)
            Xi[j-1,i] = c*τ1i + (sr*τ2i + si*τ2r)
            nτ1r = c*τ2r - (sr*τ1r + si*τ1i)
            nτ1i = c*τ2i - (sr*τ1i - si*τ1r)
            τ1r = nτ1r; τ1i = nτ1i
        end
        Xr[m,i] = τ1r; Xi[m,i] = τ1i
    end
    return
end

@inline function _solve_tturbo!(Bc, Xr, Xi, ur, ui, cs, H, Br, Bi, μ)
    copyto!(Xr, Br); copyto!(Xi, Bi)
    ldiv2_tturbo!(ur, ui, cs, H, Xr, Xi, μ)
    @inbounds @simd for j in eachindex(Bc)
        Bc[j] = complex(Xr[j], Xi[j])
    end
    return Bc
end

# ===========================================================================
# Variants. All operate on PRE-FACTORED matrices (A=F.H, C=complex(sys.C*Q),
# B=Q\sys.B, D) so the timed region isolates the per-frequency loop strategy.
# `R` is (ny, nu, nw) ComplexF64. Continuous-time frequency map assumed.
# ===========================================================================

# ---------------------------------------------------------------------------
# Modal / pole-residue formulation. This is the algorithm a Reactant.jl / XLA / GPU
# approach would actually compile: eigendecompose A = V Λ V⁻¹ ONCE, then every frequency
# is a pure elementwise broadcast + reduction (no per-frequency linear solve):
#     R(ω) = D - Σⱼ resⱼ/(λⱼ + iω),   resⱼ = (C·V)ⱼ · (V⁻¹·B)ⱼ
# Faster but less numerically robust than the Hessenberg solve for non-normal/defective A
# (accuracy depends on cond(V)). Inputs: λ (nx complex eigenvalues), res (nx residues), D.
# ---------------------------------------------------------------------------
function freqresp_modal!(R, λ, res, D, w)
    nx = length(λ)
    d = D[1, 1]
    @inbounds for i in eachindex(w)
        s = complex(0.0, w[i])
        acc = zero(eltype(R))
        @simd for j in 1:nx
            acc += res[j] / (s - λ[j])
        end
        R[1, 1, i] = d + acc
    end
    R
end

# Fully-broadcast modal form (allocates the nx×nfreq grid) — closest to what a Reactant
# trace produces; representative of the XLA-friendly array program.
function freqresp_modal_broadcast!(R, λ, res, D, w)
    s = reshape(complex.(0.0, w), 1, :)          # 1 × nfreq
    grid = res ./ (s .- λ)                        # nx × nfreq
    R[1, 1, :] .= D[1, 1] .+ vec(sum(grid; dims = 1))
    R
end

# Threaded modal loop
function freqresp_modal_threaded!(R, λ, res, D, w)
    nx = length(λ)
    d = D[1, 1]
    Threads.@threads for i in eachindex(w)
        s = complex(0.0, w[i])
        acc = zero(eltype(R))
        @inbounds @simd for j in 1:nx
            acc += res[j] / (s - λ[j])
        end
        @inbounds R[1, 1, i] = d + acc
    end
    R
end

# V1: serial, hand-rolled transcription of the package loop (sanity vs V0 baseline)
function freqresp_serial!(R, A, C, B, D, w)
    nx = size(A, 1)
    T = eltype(R)
    Bc = Vector{T}(undef, nx)
    u  = Vector{T}(undef, nx)
    cs = Vector{Tuple{real(T),T}}(undef, nx)
    @inbounds for i in eachindex(w)
        Ri = @view R[:, :, i]
        copyto!(Ri, D)
        isinf(w[i]) && continue
        copyto!(Bc, B)
        ldiv2_orig!(u, cs, A, Bc; shift = -_freqc(w[i]))
        mul!(Ri, C, Bc, -1, 1)
    end
    R
end

# V2: threaded outer loop with Threads.@spawn over contiguous chunks; per-chunk buffers
function freqresp_threaded!(R, A, C, B, D, w)
    nx = size(A, 1)
    T = eltype(R)
    nch = max(1, Threads.nthreads())
    chunks = Iterators.partition(eachindex(w), cld(length(w), nch))
    @sync for ch in chunks
        Threads.@spawn begin
            Bc = Vector{T}(undef, nx)
            u  = Vector{T}(undef, nx)
            cs = Vector{Tuple{real(T),T}}(undef, nx)
            @inbounds for i in ch
                Ri = @view R[:, :, i]
                copyto!(Ri, D)
                isinf(w[i]) && continue
                copyto!(Bc, B)
                ldiv2_orig!(u, cs, A, Bc; shift = -_freqc(w[i]))
                mul!(Ri, C, Bc, -1, 1)
            end
        end
    end
    R
end

# V2b: Polyester.@batch threading (low spawn overhead) over contiguous chunks
function freqresp_polyester!(R, A, C, B, D, w)
    nx = size(A, 1)
    T = eltype(R)
    nch = max(1, Threads.nthreads())
    chunks = collect(Iterators.partition(eachindex(w), cld(length(w), nch)))
    Polyester.@batch for ci in eachindex(chunks)
        ch = chunks[ci]
        Bc = Vector{T}(undef, nx)
        u  = Vector{T}(undef, nx)
        cs = Vector{Tuple{real(T),T}}(undef, nx)
        @inbounds for i in ch
            Ri = @view R[:, :, i]
            copyto!(Ri, D)
            isinf(w[i]) && continue
            copyto!(Bc, B)
            ldiv2_orig!(u, cs, A, Bc; shift = -_freqc(w[i]))
            mul!(Ri, C, Bc, -1, 1)
        end
    end
    R
end

# V2c: OhMyThreads with task-local buffers
function freqresp_omt!(R, A, C, B, D, w)
    nx = size(A, 1)
    T = eltype(R)
    tls_Bc = OhMyThreads.TaskLocalValue{Vector{T}}(() -> Vector{T}(undef, nx))
    tls_u  = OhMyThreads.TaskLocalValue{Vector{T}}(() -> Vector{T}(undef, nx))
    tls_cs = OhMyThreads.TaskLocalValue{Vector{Tuple{real(T),T}}}(() -> Vector{Tuple{real(T),T}}(undef, nx))
    OhMyThreads.@tasks for i in eachindex(w)
        Bc = tls_Bc[]; u = tls_u[]; cs = tls_cs[]
        Ri = @view R[:, :, i]
        copyto!(Ri, D)
        if !isinf(w[i])
            copyto!(Bc, B)
            ldiv2_orig!(u, cs, A, Bc; shift = -_freqc(w[i]))
            mul!(Ri, C, Bc, -1, 1)
        end
    end
    R
end

# Helper: solve one frequency with the split @turbo ldiv2, writing complex result into Bc
@inline function _solve_turbo!(Bc, Xr, Xi, ur, ui, cs, H, Br, Bi, μ)
    copyto!(Xr, Br); copyto!(Xi, Bi)
    ldiv2_turbo!(ur, ui, cs, H, Xr, Xi, μ)
    @inbounds @simd for j in eachindex(Bc)
        Bc[j] = complex(Xr[j], Xi[j])
    end
    return Bc
end

# V3: serial loop using the @turbo (real/imag-split) ldiv2
function freqresp_turbo!(R, A, C, B, D, w)
    nx = size(A, 1)
    T  = eltype(R)
    H  = A.data
    Br = real.(B); Bi = imag.(B)
    Bc = Vector{T}(undef, nx)
    Xr = Matrix{Float64}(undef, nx, 1); Xi = Matrix{Float64}(undef, nx, 1)
    ur = Vector{Float64}(undef, nx);    ui = Vector{Float64}(undef, nx)
    cs = Vector{Tuple{Float64,ComplexF64}}(undef, nx)
    @inbounds for i in eachindex(w)
        Ri = @view R[:, :, i]
        copyto!(Ri, D)
        isinf(w[i]) && continue
        _solve_turbo!(Bc, Xr, Xi, ur, ui, cs, H, Br, Bi, -_freqc(w[i]))
        mul!(Ri, C, Bc, -1, 1)
    end
    R
end

# V3t: serial outer loop, @tturbo (intra-loop threaded) ldiv2
function freqresp_tturbo!(R, A, C, B, D, w)
    nx = size(A, 1)
    T  = eltype(R)
    H  = A.data
    Br = real.(B); Bi = imag.(B)
    Bc = Vector{T}(undef, nx)
    Xr = Matrix{Float64}(undef, nx, 1); Xi = Matrix{Float64}(undef, nx, 1)
    ur = Vector{Float64}(undef, nx);    ui = Vector{Float64}(undef, nx)
    cs = Vector{Tuple{Float64,ComplexF64}}(undef, nx)
    @inbounds for i in eachindex(w)
        Ri = @view R[:, :, i]
        copyto!(Ri, D)
        isinf(w[i]) && continue
        _solve_tturbo!(Bc, Xr, Xi, ur, ui, cs, H, Br, Bi, -_freqc(w[i]))
        mul!(Ri, C, Bc, -1, 1)
    end
    R
end

# V4b: Polyester.@batch + @turbo ldiv2 (combination of the two best individual strategies)
function freqresp_polyester_turbo!(R, A, C, B, D, w)
    nx = size(A, 1)
    T  = eltype(R)
    H  = A.data
    Br = real.(B); Bi = imag.(B)
    nch = max(1, Threads.nthreads())
    chunks = collect(Iterators.partition(eachindex(w), cld(length(w), nch)))
    Polyester.@batch for ci in eachindex(chunks)
        ch = chunks[ci]
        Bc = Vector{T}(undef, nx)
        Xr = Matrix{Float64}(undef, nx, 1); Xi = Matrix{Float64}(undef, nx, 1)
        ur = Vector{Float64}(undef, nx);    ui = Vector{Float64}(undef, nx)
        cs = Vector{Tuple{Float64,ComplexF64}}(undef, nx)
        @inbounds for i in ch
            Ri = @view R[:, :, i]
            copyto!(Ri, D)
            isinf(w[i]) && continue
            _solve_turbo!(Bc, Xr, Xi, ur, ui, cs, H, Br, Bi, -_freqc(w[i]))
            mul!(Ri, C, Bc, -1, 1)
        end
    end
    R
end

# V4: threaded (@spawn) + @turbo ldiv2
function freqresp_threaded_turbo!(R, A, C, B, D, w)
    nx = size(A, 1)
    T  = eltype(R)
    H  = A.data
    Br = real.(B); Bi = imag.(B)
    nch = max(1, Threads.nthreads())
    chunks = Iterators.partition(eachindex(w), cld(length(w), nch))
    @sync for ch in chunks
        Threads.@spawn begin
            Bc = Vector{T}(undef, nx)
            Xr = Matrix{Float64}(undef, nx, 1); Xi = Matrix{Float64}(undef, nx, 1)
            ur = Vector{Float64}(undef, nx);    ui = Vector{Float64}(undef, nx)
            cs = Vector{Tuple{Float64,ComplexF64}}(undef, nx)
            @inbounds for i in ch
                Ri = @view R[:, :, i]
                copyto!(Ri, D)
                isinf(w[i]) && continue
                _solve_turbo!(Bc, Xr, Xi, ur, ui, cs, H, Br, Bi, -_freqc(w[i]))
                mul!(Ri, C, Bc, -1, 1)
            end
        end
    end
    R
end
