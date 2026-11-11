# Diagram compiler (`fepdmc_toolkit.diagram_compiler`, experimental)

A compiler for fixed-order electron–phonon diagram families. It makes the
diagram sum explicit as a serializable `DiagramIR`, lowers it to a shared DAG
and an exact evaluator, derives a cheaper evaluator from the same IR, pairs the
two in an exact delayed-acceptance kernel, lowers the kernel to a typed
control-flow IR (`NativeKernelIR`) and emits Fortran from it.

```text
R1 PlanSpec ──import──> DiagramIR ──> exact DAG ──> exact evaluator ─┐
                            │                                         ├─> delayed-acceptance kernel
                            └───────> cheap DAG ──> cheap evaluator ──┘            │
                                                                       NativeKernelIR ──> Fortran
```

Supported families: scalar Holstein (validated orders n = 1–4) and Hermitian
two-band (n = 1–3). The physical object is `F_n(x)`, the sum over all
`(2n−1)!!` perfect chord matchings of `2n` vertices at one sample `x`.

## Quick start

```python
from fepdmc_toolkit.diagram_compiler import build_scalar_ir, compile_evaluator, compile_cheap
from fepdmc_toolkit.diagram_compiler.evaluator import Binding, torus_point

ir = build_scalar_ir(2)  # scalar Holstein, order 2
L = 4
x = Binding(
    k=torus_point(0, L),
    q=(torus_point(1, L), torus_point(2, L)),
    tau=(0.05, 0.12, 0.31, 0.47),
    t=1.0,
    omega=0.8,
    g=0.5,
    L=L,
)
print(compile_evaluator(ir)(x)["F"])  # exact diagram sum
print(compile_cheap(ir).evaluate(x)["F_hat"])  # cheap evaluator
```

`examples/diagram_compiler/quickstart.py` runs a complete example. Each
stage below has a script, `scripts/dc_stage<N>_….py`, that repeats its checks,
writes the records and exits with status 1 if a check fails.

## Stage 1: exact compiler

**Physics.** Scalar: `G(p, Δτ) = exp(−2t(1 − cos p) Δτ)`,
`D(a, b) = exp(−ω (τ_b − τ_a))`,
`F_n = (g/√L)^(2n) Σ_C Π G Π D`; every diagram has sign +1 and multiplicity 1.
Two-band: `H(k) = [[ξ_k, δ], [δ, ξ_k + gap]]`, `U = exp(−H Δτ)`,
`M = (g/√L) diag(1, −1)` at every vertex; traversing the vertices forward,
`state = M @ state`, then `state = U @ state`; the final matrix times the
phonon factors summed over diagrams gives `Op`, and `F = trace(Op)`. Chords
receive q slots by first-opening rank; an internal electron segment carries
`k − Σ q[open slots]`. No external legs, Jacobians, factorials, tadpoles or
dressed propagators are added.

**IR.** `DiagramIR` (schema 1.0.0) holds the family (model, order, bands,
dispersion, phonon and normalization conventions, rule version, global
factors), the variables (external k, q and time slots, strict time order,
dynamic parameter names) and the canonical diagram list with explicit
factors (electron interval and momentum form, phonon endpoints, vertex,
prefactor). Dynamic values are `(k, q, τ, t, ω, g, L)` plus `(δ, gap)` for two
bands; order, topology, signs and multiplicities are structural. JSON field
order is deterministic, and deserialization rejects foreign schemas, unknown
fields, non-integer indices, duplicate or missing diagrams and wrong factors.
Compilation snapshots the IR.

**Canonicalization and sharing.** Vertices keep physical time order (no
sorting on the caller's behalf); chords are earlier-to-later; electron keys
contain the interval, open-slot set and integer momentum form; matrix
products are never reordered. The DAG is built by ordered structural
hash-consing: identical leaves and factor-chain prefixes are shared, and a
node with several uses is evaluated once per binding. There is no
floating-point merging or sample-dependent folding.

| family / order | diagrams | naive nodes | unique DAG nodes | arithmetic ops | R1 nodes |
|---|---:|---:|---:|---:|---:|
| scalar 1 | 1 | 7 | 7 | 3 | 7 |
| scalar 2 | 3 | 35 | 28 | 15 | 28 |
| scalar 3 | 15 | 257 | 124 | 91 | 95 |
| scalar 4 | 105 | 2417 | 803 | 732 | 315 |
| two-band 1 | 1 | 12 | 12 | 6 | 11 |
| two-band 2 | 3 | 62 | 47 | 27 | 43 |
| two-band 3 | 15 | 452 | 194 | 142 | 142 |

The DAG does not recover all of the recursive R1 evaluator's distributive
sharing (compare the last two columns); node counts are not timings.
`lower_grouped` is a separate hand-written transcription of the R1
recurrence, used for comparison only.

**Validation.** Tolerance `|value − reference| ≤ 1e-12 + 1e-10 |reference|`,
non-finite results rejected. References: the R1 recursive builder
(`fepdmc_toolkit.scheme_d`, executed by its own interpreter) and an
independent explicit enumeration (`index_oracle`). 224 fixed bindings over
all supported families and orders agree with both; the maximum absolute
difference from R1, including every two-band `Op` entry, is 2.78 × 10⁻¹⁷.
Seeded randomized tests vary L, k, q, times and parameters; asymmetric
controls detect sign, momentum, phonon-endpoint, sharing, q-slot and
non-commutative multiplication faults; isolated subprocesses evaluate JSON IR
with the R1 code blocked.

## Stage 2: cheap evaluator

Policy `propagator_only_v1` (`validation/diagram_compiler/CHEAP_POLICY.json`):

| action | categories |
|---|---|
| keep | electron and phonon propagators, prefactor, sign, multiplicity, global factors, all pairings, q terms |
| drop | vertices, `matmul`, `trace` |
| approximate | two-band matrix propagator → scalar `Gel` of the same momentum; empty two-band global factors → `(g/√L)^(2n)` |
| refuse | any other category |

The cheap evaluator is a state weight `W_hat(x)`; the transition score is
`ell_hat(x, y) = log W_hat(y) − log W_hat(x)`, antisymmetric by
construction. `compile_cheap` emits a scalar DAG (`input`, `const`, `add`,
`mul`) and never calls the exact lowering or the R1 code. On two bands the
cheap DAG has no vertex, eigen-decomposition, `matmul` or `trace` (7, 28, 124
unique nodes against 12, 47, 194 exact). On the scalar family the policy
equals the exact evaluator, so it is not cheaper there.

```python
exact = compile_exact(ir)
cheap = compile_cheap(ir, policy="propagator_only_v1")
cheap.transition_score(x, y)
```

## Stage 3: delayed-acceptance kernel

Target `π(x) = F(x)` on the set where F is finite, real (|Im F| ≤ 1e-12) and
strictly positive (`positive_real_F_v1`); `|F|`, `F²` and complex moduli are
refused. Proposals are `symmetric` or carry a caller-supplied
`log q(y, x) − log q(x, y)`.

Christen–Fox two-stage kernel, no clipping:
`α1 = min(1, exp(ell_hat))`,
`α2 = min(1, exp(ell_R − ell_hat))` with
`ell_R = log π(y) − log π(x) + log q(y, x) − log q(x, y)`. Then
`α1 α2 / (α1 α2)_reverse = π(y) q(y, x) / (π(x) q(x, y))` exactly.
A stage-1 reject evaluates only the cheap weights; a stage-1 pass evaluates
`exact(y)` and `exact(x)` (or uses a correctly keyed cache for `exact(x)`).

```python
kernel = compile_delayed_acceptance(ir)
result = kernel.evaluate_transition(x, y, log_q_reverse_minus_forward, u1, u2)
```

On a three-state chain with a different cheap weight and asymmetric `q`, the
detailed-balance residual is 1.39 × 10⁻¹⁷; mutations (no stage 2, missing
`−ell_hat`, inverted Hastings factor, cheap used as exact) break balance.

## Stage 4: NativeKernelIR

A typed SSA control-flow IR (`bool`, `i64`, `f64`, `c128`, fixed 2×2
`c128_m2`) with a binding layout, reusable cheap and exact graphs, and blocks
with `br`, `br_cond` and `return`. Primitive set `native_kernel_v1`: electron
(scalar and two-band), phonon, prefactor, σ_z vertex, `matmul`, `trace`,
`scale_m2`, and control operations (`log_positive_real`, `log_uniform`,
`min_f64`, `lt_f64`, `check_log_q`, `call_graph`, `phi`). Status codes:
`OK = 0`, `INVALID_EXACT_TARGET = 1`, `INVALID_CHEAP_WEIGHT = 2`,
`INVALID_PROPOSAL_RATIO = 3`, `INVALID_RANDOM_UNIFORM = 4`,
`CACHE_KEY_MISMATCH = 5`, `INVALID_IR = 6`. No exact operation is reachable
on the stage-1 reject path. `interpret_native_kernel` executes the IR without
the Python kernel and matches stage 3 on all paths. The two-band n = 1 IR is
recorded in `validation/diagram_compiler/NATIVE_KERNEL_IR.json`.

## Stage 5: Fortran backend

`fortran_codegen.py` emits a kernel-specific subroutine from NativeKernelIR
only; `fortran/native_kernel_runtime_v1.f90` implements the primitives
(two-band propagator via `cosh`/`sinh`, no LAPACK). Types map to
`iso_fortran_env` kinds; the CFG becomes labeled blocks; kernels set `status`
and return instead of stopping. ABI for two-band n = 1:

```text
subroutine evaluate_da_kernel( &
    x_k, x_q(1), x_tau(2), x_t, x_omega, x_g, x_L, x_delta, x_gap, &
    y_k, y_q(1), y_tau(2), y_t, y_omega, y_g, y_L, y_delta, y_gap, &
    log_q_reverse_minus_forward, u1, u2, &
    exact_x_valid, exact_log_weight_x, &
    accepted, reject_stage, ell_hat, ell_R_valid, ell_R, &
    exact_y_evaluated, status)
```

The caller owns the random numbers and sets `exact_x_valid` only when the
cached log weight belongs to `x`. The same IR always yields the same source:
`validation/diagram_compiler/generated_twoband_n1.f90`, SHA256 `ff4aa2de…`.
Compiled kernels match the interpreter for every path, status and primitive,
and on a stage-1 reject no exact primitive is called.

## Stage 6: native toolchain validation

The generated two-band n = 1 kernel and an emitted driver
(`dc_validate.f90`) are compiled with gfortran 11.4 on Ubuntu 22.04
(`-O0 -std=f2008 -ffree-line-length-none -fcheck=bounds`) as a separate
executable next to an unmodified clone of the pinned FEP-DMC tree, which is
checked to be pristine. Fixtures come from the IR interpreter. All cases
match: stage-1 reject, stage-2 reject, accept, asymmetric `log q`, cached and
uncached `exact(x)`, and the error statuses for `g = 0`, non-finite `log q`
and invalid uniforms; discrete fields exactly and floating-point fields with
zero error (`validation/diagram_compiler/stage6_native_differential.csv`).
Production FEP-DMC does not call the generated kernel.

To run it:

```bash
docker build -f integration/diagram_compiler/Dockerfile.u22 -t fepdmc-toolkit-u22:local .
git clone https://github.com/yaoluo/FEP-DMC /tmp/fep-dmc
git -C /tmp/fep-dmc checkout 05d08449cffdbd0dfbbbf5009add5cc887bc754b
python scripts/dc_stage6_fepdmc_toolchain.py --fep-dmc /tmp/fep-dmc --image fepdmc-toolkit-u22:local
```

## Relation to P1

Both use the same two-stage algebra, but the objects differ: the compiler
kernel targets the positive diagram sum `F` with a state-weight cheap score
and an explicit proposal ratio; P1 targets the native `abs(Re M) · P_kchange`
of the swap update with a clipped local transition score. The compiler does
not replace P1.

## Limits

Two hard-coded families; orders up to 6 are admitted as a resource bound,
but only the orders above are validated. Pairing enumeration is factorial.
No variable-order compilation, topology sampler, general material provider,
GPU or LLVM backend, and no speed-up claim. A new family needs its own
external-leg, band, gauge, sign and multiplicity rules and its own
independent oracle before it can be admitted.
