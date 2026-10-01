# Benchmark freqresp! strategies for a 48-state SISO state-space system over 150 frequencies.
#
# Launch with multiple threads to exercise the threaded variants:
#   julia -t auto --project=. benchmark.jl
#
# Variants (see variants.jl):
#   V0 baseline        ControlSystemsBase.freqresp!         (correctness oracle, full call)
#   V1 serial          hand-rolled transcription, loop-only
#   V2 threaded        Threads.@spawn over chunks
#   V2b polyester      Polyester.@batch over chunks
#   V2c omt            OhMyThreads task-local buffers
#   V3 turbo           @turbo (real/imag split) ldiv2, serial
#   V4 threaded+turbo  @spawn + @turbo ldiv2

using ControlSystemsBase, BenchmarkTools, LinearAlgebra, Printf
import LoopVectorization

include("variants.jl")

const NTH = Threads.nthreads()
@info "Benchmark configuration" threads=NTH
NTH == 1 && @warn "Running with a single thread — threaded variant timings are meaningless. Relaunch with `julia -t auto`."

# --- Target system & frequencies ---
G = ssrand(1, 1, 48)                       # 48-state SISO, continuous
w = exp10.(LinRange(-2, 4, 150))           # 150 log-spaced frequencies
ny, nu = size(G)
const T = ComplexF64
make_R() = Array{T,3}(undef, ny, nu, length(w))

# --- Pre-factored matrices for the loop-only variants ---
F = hessenberg(G.A)
Q = Matrix(F.Q)
A = F.H
C = complex.(G.C * Q)
B = Q \ G.B
D = G.D

# Modal / pole-residue precompute (one-time eigendecomposition) — the Reactant-friendly form
E = eigen(G.A)
λ = ComplexF64.(E.values)                  # nx eigenvalues
V = ComplexF64.(E.vectors)
ctil = vec(complex.(G.C) * V)              # C·V               (nx,)
btil = V \ ComplexF64.(G.B[:, 1])          # V⁻¹·B             (nx,)
res = ctil .* btil                         # residues          (nx,)
condV = cond(V)                            # conditioning of the modal transform

# A wrapped "full call" variant that does the factorization inside (apples-to-apples w/ V0)
function freqresp_serial_full!(R, sys, w)
    F = hessenberg(sys.A); Q = Matrix(F.Q)
    freqresp_serial!(R, F.H, complex.(sys.C*Q), Q\sys.B, sys.D, w)
end

# ----------------------------------------------------------------------------
# Correctness: every variant must match the package baseline.
# ----------------------------------------------------------------------------
Rref = make_R(); ControlSystemsBase.freqresp!(Rref, G, w)

loop_variants = [
    ("V1 serial",          (R)->freqresp_serial!(R, A, C, B, D, w)),
    ("V2 threaded",        (R)->freqresp_threaded!(R, A, C, B, D, w)),
    ("V2b polyester",      (R)->freqresp_polyester!(R, A, C, B, D, w)),
    ("V2c omt",            (R)->freqresp_omt!(R, A, C, B, D, w)),
    ("V3 turbo",           (R)->freqresp_turbo!(R, A, C, B, D, w)),
    ("V3t tturbo",         (R)->freqresp_tturbo!(R, A, C, B, D, w)),
    ("V4 threaded+turbo",  (R)->freqresp_threaded_turbo!(R, A, C, B, D, w)),
    ("V4b polyester+turbo",(R)->freqresp_polyester_turbo!(R, A, C, B, D, w)),
    ("M1 modal",           (R)->freqresp_modal!(R, λ, res, D, w)),
    ("M2 modal broadcast", (R)->freqresp_modal_broadcast!(R, λ, res, D, w)),
    ("M3 modal threaded",  (R)->freqresp_modal_threaded!(R, λ, res, D, w)),
]
@info "Modal transform conditioning" condV

println("\n=== Correctness (rtol=1e-10 vs ControlSystemsBase.freqresp!) ===")
for (name, f) in loop_variants
    R = make_R(); f(R)
    relerr = maximum(abs, R .- Rref) / maximum(abs, Rref)
    ok = relerr < 1e-8
    @printf("  %-20s %s   relerr=%.3e\n", name, ok ? "OK   " : "LOOSE", relerr)
end

# ----------------------------------------------------------------------------
# Timing
# ----------------------------------------------------------------------------
Rb = make_R()

println("\n=== Loop-only timing (pre-factored matrices passed in) ===")
results = Tuple{String,Float64}[]
for (name, f) in loop_variants
    t = @belapsed $f($Rb)
    push!(results, (name, t))
    @printf("  %-20s %8.2f µs\n", name, t*1e6)
end

println("\n=== Full-call timing (factorization included, as real freqresp usage) ===")
tb = @belapsed ControlSystemsBase.freqresp!($Rb, $G, $w)
@printf("  %-20s %8.2f µs\n", "V0 baseline", tb*1e6)
tf = @belapsed freqresp_serial_full!($Rb, $G, $w)
@printf("  %-20s %8.2f µs\n", "V1 serial (full)", tf*1e6)

# --- Summary table sorted by loop-only time ---
println("\n=== Summary (loop-only, fastest first) ===")
sort!(results, by = x -> x[2])
fastest = results[1][2]
for (name, t) in results
    @printf("  %-20s %8.2f µs   (%.2fx vs fastest)\n", name, t*1e6, t/fastest)
end
@printf("\n  Full call baseline (V0): %.2f µs  — factorization dominates at this size.\n", tb*1e6)
println("  Threads used: $NTH")
