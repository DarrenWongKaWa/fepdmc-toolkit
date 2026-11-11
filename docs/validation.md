# Validation

What each component has been checked against, and where the record is.
"Source level" means the transform is verified on the pinned FEP-DMC tree
by exact file hashes, without compiling. "Native" means a compiled
`perturbo.x` or generated Fortran was run.

## Test suite (no native code)

`pytest` runs in well under a minute without FEP-DMC, material data or
Docker (the tests that compile Fortran are skipped without gfortran). Main
groups:

| tests | checks |
|---|---|
| `test_da_balance.py`, `test_forward_reverse.py`, `test_audit_guards.py` | delayed acceptance restores the native `abs(Re M)` ratio (`|M|` is a different target); occupancy reverse is antisymmetric; invalid inputs fail |
| `test_energy_zero.py`, `test_statistics.py`, `test_c5_formation_q.py`, `test_reject_state.py` | Q = ratio − E_bare; recorded statistics, HAC failures, counters and classification |
| `test_p1_clean_source.py` | the benchmarked P1 Fortran: dump and reverse checks gated, stage-1 return before `cal_gkq`, counters independent of dumps |
| `test_fep_dmc_upstream.py` | upstream pin and the fail-closed verifier |
| `test_c0_adapter.py`, `test_p1_public_adapter.py`, `test_p1_public_fortran_contract.py` | C0 and P1 apply/verify state machines, composition in either order, refusal of unknown trees, forbidden symbols and falsifiers in the P1 module |
| `test_fepdmc_tools.py` | patch registry and anchors, refusal of unknown trees by `prepare`, CLI, run driver (dry run, checks before start, one container per chain), the pooled estimator and the exactness check |
| `test_analysis.py` | every analysis on deterministic stand-ins for run directories, against recorded reference outputs; missing input is reported |
| `test_prepare_upstream.py` | `prepare` on a clone of the pinned commit reproduces `provenance/PREPARED_TREE.json` for both profiles, writes nothing in a dry run and refuses a patched tree (skipped without the clone) |
| `test_r1_oracles.py`, `test_r1_grouped.py` | R1 metadata; the grouped toy against its dense oracle |
| `test_diagram_ir.py` … `test_native_kernel.py`, `test_fortran_*.py` | the diagram compiler: IR validation, differential tests against R1 and the index oracle, DAG sharing, cheap evaluator, delayed-acceptance kernel, NativeKernelIR, Fortran backend (compiled with gfortran) |
| `test_dc_golden.py`, `test_dc_fep_dmc_validate.py` | the generated Fortran, NativeKernelIR and runtime are the natively validated artifacts; Docker-based native check (skipped without the image) |
| `test_integrity.py` | validated artifacts match `validation/SHA256SUMS` |
| `test_packaging.py` | single version definition, sdist contents, no private paths |

The continuous integration also runs `ruff`, the tests from an unpacked
source distribution, the C0 and P1 adapters and `prepare` on a clone of the
pinned commit, and the six validation stages of the diagram compiler. Each
stage script exits with status 1 if one of its checks fails.

## FEP-DMC patches

- Every move passes a round trip C → C′ → C that restores log w and the
  proposal probabilities to ≤ 4 × 10⁻¹⁵, plus per-chain momentum, gauge and
  sign self-checks.
- `fepdmc-toolkit validate` runs three independent kernels (fixed native,
  toolkit moves only, both) on the same target at low order. At maxOrder 7
  they agree; the fixed build reproduced the total energy E = −0.8763 … −0.8767
  for LiF-hole at T = 500 K.
- Effects of the fixes and the move efficiencies were measured on LiF-hole,
  LiF-electron and SrTiO₃ ([fepdmc.md](fepdmc.md), records in
  [validation/fepdmc/](../validation/fepdmc/README.md)).

## C0 adapter

Source level: pristine pre-image and exact post-image hashes at the pinned
commit ([provenance/C0_PUBLIC_PATCH.json](../provenance/C0_PUBLIC_PATCH.json)).
The same one-line fix was used natively as `--wqfix`.

## P1 adapter

- Source level: exact hashes for the four states (pristine, C0, P1, both) and
  the falsification tests listed in [p1-adapter-design.md](p1-adapter-design.md).
- Native: built from the pinned tree and run on a short public LiF JJ
  fixture; P1 off reproduces the pristine counts; stage-1 rejects skip
  `cal_gkq*`; counters close ([p1.md](p1.md),
  [validation/p1/](../validation/p1/NATIVE_VALIDATION.json)).

## Benchmarked P1 build

Event-level fixture on LiF: 4372 stage-1 rejects with zero recorded
electron–phonon and environment increments; 6350 accepted-swap reverse checks
passed. The fingerprint used to compare live state is partial; the stronger
evidence for the skip is the recorded zero increment. Benchmark in
[results.md](results.md).

## Diagram compiler

- Stage 1: 224 fixed bindings (scalar Holstein n = 1–4, Hermitian two-band
  n = 1–3) agree with the R1 reference and an independent index oracle; max
  absolute error including every two-band `Op` entry 2.78 × 10⁻¹⁷ at
  `atol = 1e-12`, `rtol = 1e-10`
  ([stage1_differential.csv](../validation/diagram_compiler/stage1_differential.csv)).
- Stage 3: detailed-balance residual 1.39 × 10⁻¹⁷ on a three-state chain.
- Stage 5: the compiled kernel matches the NativeKernelIR interpreter for all
  paths and statuses; every primitive matches its oracle.
- Stage 6: the generated two-band n = 1 kernel, compiled with gfortran 11.4
  on Ubuntu 22.04 next to the unmodified pinned FEP-DMC tree, matches the
  interpreter on all fixtures; discrete fields exactly, floating-point fields
  with zero error
  ([stage6_native_differential.csv](../validation/diagram_compiler/stage6_native_differential.csv),
  recorded on aarch64; the check also passes on x86_64).

Details in [diagram-compiler.md](diagram-compiler.md).

## Not validated

See [limitations.md](limitations.md).
