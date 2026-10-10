# `freqresp!` optimization benchmark

Standalone benchmark of strategies for `freqresp!` on the target use case:
a **48-state SISO continuous state-space system over 150 frequencies**.

Nothing in the package source is modified — all candidate implementations live in
`variants.jl`. The package's own `freqresp!` is the correctness oracle.

## Run

```bash
julia --project=. -e 'using Pkg; Pkg.develop(path="../lib/ControlSystemsBase"); \
  Pkg.add(["BenchmarkTools","LoopVectorization","Polyester","OhMyThreads"]); \
  Pkg.instantiate(); Pkg.precompile()'

JULIA_EXCLUSIVE=1 julia -t auto --project=. benchmark.jl   # threads fixed at startup
```

## Variants

| ID  | Strategy                                                              |
|-----|----------------------------------------------------------------------|
| V0  | `ControlSystemsBase.freqresp!` (oracle, full call incl. factorization) |
| V1  | serial hand-rolled loop, original `ldiv2!` (loop-only)               |
| V2  | `Threads.@spawn` over contiguous frequency chunks                    |
| V2b | `Polyester.@batch` over chunks                                       |
| V2c | `OhMyThreads.@tasks` with `TaskLocalValue` buffers                   |
| V3  | serial loop, `@turbo` `ldiv2!` (real/imag split)                     |
| V4  | `@spawn` + `@turbo` `ldiv2!`                                          |
| V4b | `Polyester.@batch` + `@turbo` `ldiv2!`                               |

All variants operate on pre-factored matrices (`A=F.H`, `C=complex.(sys.C*Q)`,
`B=Q\sys.B`, `D`) so the "loop-only" timing isolates the per-frequency strategy from the
one-time Hessenberg factorization. All pass the `isapprox(·, ·; rtol=1e-10)` correctness
gate against the package baseline.

## The `@turbo` / complex caveat

LoopVectorization's `@turbo` does **not** support `Complex`. A bare `@turbo` on the
complex `ldiv2!` loops will not compile, and a naive `reinterpret`-to-real is *silently
wrong* (complex multiply has re/im cross terms). `ldiv2_turbo!` therefore carries the
state as separate real/imag `Float64` arrays and writes out the complex arithmetic by
hand, so the two hot `j = 1:k-2` loops (X update, u update) become valid real `@turbo`
loops. Result matches the baseline to ~1e-14.

## Results (this machine: 24 logical cores)

Loop-only time, µs (lower is better):

| Variant            | -t 4  | -t 8  | -t 24 |
|--------------------|-------|-------|-------|
| V1 serial          | 545   | 566   | 565   |
| V3 turbo           | 360   | 372   | 372   |
| V2 threaded        | 147   |  90   | 106   |
| V2b polyester      | 203   | 108   |  95   |
| V2c omt            | 154   |  96   | 127   |
| V4 threaded+turbo  | **104** | **78** | 101 |
| V4b polyester+turbo| 142   |  80   | **93** |

Full-call baseline (V0, factorization included): ~625–655 µs.

## Findings

1. **The per-frequency loop, not the factorization, dominates.** The serial loop is
   ~565 µs while the one-time Hessenberg factorization is only ~85 µs (full call ≈ loop +
   factorization ≈ 650 µs). So optimizing the loop is worthwhile.

2. **`@turbo` is a solid, free win: ~1.5×** (565 → ~370 µs serial), more than expected for
   such short loops, with no threads and negligible accuracy loss. It does require the
   manual real/imag split.

3. **Threading gives ~6–7×** but scales sub-linearly: at this size **8 threads beats 24**
   (150 freqs / 24 threads ≈ 6 freqs each — task overhead and oversubscription dominate).
   The sweet spot here is ~8 threads.

4. **Best overall is threading + `@turbo`** (~78 µs at 8 threads, a **~7.2×** speedup over
   serial). `@spawn`+`@turbo` (V4) is consistently near-best across thread counts;
   `Polyester`+`@turbo` (V4b) edges ahead only at 24 threads.

5. **Scheduler choice is secondary.** `@spawn`, `Polyester`, and `OhMyThreads` are within
   ~1.4× of each other; `@spawn` is the most robust across thread counts and adds no
   dependency.

## Reactant.jl and the modal reformulation

Reactant.jl traces Julia into MLIR/XLA and wants vectorized array programs. The custom
Givens-based `ldiv2!` (scalar indexing, data-dependent control flow, `givensAlgorithm`)
does **not** trace into efficient XLA. To use Reactant you must first rewrite the
computation as batched array ops — and the natural rewrite is the **modal / pole-residue**
form: eigendecompose `A = V Λ V⁻¹` once, then

```
R(ω) = D + Σⱼ resⱼ / (iω − λⱼ),   resⱼ = (C·V)ⱼ · (V⁻¹·B)ⱼ
```

Every frequency is then a pure broadcast + reduction over an `nx × nfreq` grid — no
per-frequency linear solve. Benchmarked in plain Julia (variants M1/M2/M3):

| Variant            | µs (8 threads) |
|--------------------|----------------|
| M1 modal serial    | 42             |
| M2 modal broadcast | 59             |
| M3 modal threaded  | **10**         |

So the modal form alone is **13×** (serial) to **54×** (threaded) faster than the current
serial Hessenberg loop, and beats the best Hessenberg micro-optimization (threaded+`@turbo`,
~73 µs) by 7×. Accuracy matched the baseline to ~1e-14 here (`cond(V)=147`).

**Verdict on Reactant for this problem:** not worth it. The data is tiny (48×150 ≈ 7200
complex numbers) and the modal form already runs in ~10 µs on CPU — below the per-call XLA
dispatch / GPU kernel-launch latency (tens of µs). Reactant would add a heavy dependency
and long compile times for no gain at this size. It would only pay off at much larger scale
(thousands of states, very many frequencies, or batches of systems) or on GPU — and there
the same modal array program is what you would feed it.

### When the modal form degrades or fails

The modal form's accuracy is governed by `cond(V)` (the conditioning of the eigenvector
matrix). Checked against a 256-bit `BigFloat` resolvent reference:

| System                              | cond(V) | Hessenberg relerr | modal relerr |
|-------------------------------------|---------|-------------------|--------------|
| random (near-normal), nx=20         | 2e1     | 2e-15             | 8e-15 ✓      |
| 8 repeated poles at −1              | 1.5e110 | 1e-15             | **1.0** ✗    |
| 12 repeated poles at −1             | 7.6e172 | 2e-15             | **1.0** ✗    |
| clustered eigenvalues (non-normal)  | 1e20    | 8e-16             | **8.7e2** ✗  |
| lightly damped, eval across pole    | 1e0     | 6e-15             | 2e-14 ✓      |

Failure conditions, in order of severity:

1. **Defective (non-diagonalizable) `A`** — a repeated eigenvalue with a non-trivial Jordan
   block. `A = VΛV⁻¹` does not exist (`V` singular), and the simple-pole form cannot
   represent the response (needs `rⱼ/(s−λ)ᵏ` terms). Total failure (relerr ≈ 1).
2. **Near-defective / strongly non-normal `A`** — diagonalizable but `cond(V)` huge.
   Forming `V⁻¹B` and the residue sum loses ≈ `log₁₀(cond(V))` digits with catastrophic
   cancellation; error scales with `cond(V)`.
3. **Clustered eigenvalues** — nearly-parallel eigenvectors → ill-conditioned `V` (a milder
   case 2).

These are common in control (repeated integrators, cascaded identical lags,
Butterworth/Bessel clustered poles). Note that evaluating *near a pole* (`iω ≈ λⱼ`) is **not**
a modal-specific weakness — the response genuinely blows up there for every method; with a
well-conditioned `V` the modal form is as accurate as Hessenberg across the resonance.

**Why the Hessenberg solve is immune:** it computes `(iωI − A)⁻¹B` via a backward-stable
linear solve whose accuracy depends on `cond(iωI − A)` (benign except genuinely at a pole),
**not** on `cond(V)`, and it never diagonalizes so it handles defective `A` fine. That
robustness is why the package defaults to it.

**Caveat / rule of thumb:** the modal fast path is safe only when `A` is near-normal with
simple, well-separated poles. Guard it by checking `cond(V)` — beyond ~`1/√eps ≈ 1e8` half
the digits are gone, and it should fall back to the Hessenberg solve. It is an excellent
opt-in fast path for well-conditioned systems, not a safe drop-in default.

## Recommendation

For the package, the lowest-risk improvement is dropping in the real/imag-split `@turbo`
`ldiv2!` (≈1.5× for free, no new threading semantics). If the loop is a known hotspot,
add an opt-in threaded path (`@spawn` over chunks with per-task buffers) — combined with
`@turbo` it reaches ~7× — but keep the default serial, since threading helps only for
many frequencies and benefits from a bounded thread count.
