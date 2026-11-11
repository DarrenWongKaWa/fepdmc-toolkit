# Names used in this repository

Short names from the development of the code remain in file names, switches
and recorded results. This page says what each one means.

## Upstream

| name | meaning |
|---|---|
| FEP-DMC | first-principles electron–phonon diagrammatic Monte Carlo, the code of Luo, Park and Bernardi built on Perturbo and Quantum ESPRESSO 6.5 |
| native | the unmodified upstream code path; with a patch applied, the path taken when no toolkit switch is set |
| pin, pinned commit | FEP-DMC commit `05d08449cffdbd0dfbbbf5009add5cc887bc754b`, the only upstream revision the patches accept ([upstream.md](upstream.md)) |
| EZ | the `diagmc-EZ` calculation mode of FEP-DMC; its energy estimator is a ratio of block sums |
| JJ | the `diagmc-JJ` calculation mode and the files named after it. `diagMC_JJ_updates.f90` holds the update routines that the patches and adapters edit (`update_swap`, `add_external_ph`, `remove_external_ph`) |
| Method 0 | `DMC_Method = 0`, the sampling method of all runs recorded here |
| external pair | a phonon line that wraps the imaginary-time boundary; added and removed by `add_external_ph` and `remove_external_ph` |

## Quantities

| name | meaning |
|---|---|
| E, `Etrue` | total energy estimate printed by the native code |
| E_bare | bare band energy; `E_BARE_LIF_METHOD0` for the LiF-electron benchmark |
| Q | polaron formation energy, Q = E − E_bare ([statistics.md](statistics.md)) |
| ⟨s⟩, `<g/Z>` | average sign of the sampled weight |
| τ_max, `tau_max` | projection length 1/kT in eV⁻¹ (232.09 at 50 K) |
| K | cost at equal precision, wall time × jackknife variance of Q |
| HAC | heteroskedasticity- and autocorrelation-consistent standard error (Bartlett kernel) |
| RB | Rao–Blackwellization: averaging an estimator over a group of configurations of the same chain |

## Patches and adapters

| name | meaning |
|---|---|
| `FUTUREB_…` | prefix of every run-time switch of the patched binary. "Future B" was the working name of the project. The name is part of the patched Fortran, which is kept byte-identical to the validated builds, so the switches keep it |
| wqfix, extfix, extrmfix | the three opt-in fixes for the external-phonon updates ([fepdmc.md](fepdmc.md)) |
| profile `fixes` | the build fixes plus the three opt-in fixes |
| profile `research` | `fixes` plus the Monte Carlo moves and the diagnostics |
| chain A | a chain of the native kernel |
| chain B, bchain | the grouped-measure chain: phonon modes summed over a non-crossing set of lines, sampled with the native update as proposal |
| chain C | in `fepdmc-toolkit compare`, the modified chain compared with chain A |
| C0 | the source adapter that inserts the phonon-frequency refresh into `add_external_ph` ([c0-adapter.md](c0-adapter.md)) |
| P1 | two-stage delayed acceptance in the swap update ([p1.md](p1.md)); `ell_hat` is the score of the cheap first stage, `ell_R` the logarithm of the exact acceptance ratio |
| `LINEAR_DA` | the run-time switch of P1 |

## Benchmark records

| name | meaning |
|---|---|
| B0, B-best | two native builds without delayed acceptance. P1 is compared with B-best. The builds are identified by the SHA256 of their binaries in `benchmarks/` |
| `benchmarks/c5_dev` | formation energy of the three builds, six chains each |
| `benchmarks/p1_vs_bbest` | timing of P1 against B-best |
| `benchmarks/c0_c1` | acceptance check of the frequency refresh and the profile of the swap update |
| `benchmarks/r1_fixed_order` | fixed-order R1 evaluator against common-subexpression elimination |

## Diagram compiler

| name | meaning |
|---|---|
| F_n(x) | sum over all (2n − 1)!! pairings of 2n vertices at one sample point x, at fixed vertex times |
| R1 | the reference evaluator of F_n: a compile-once recursion over open chords. It lives in `fepdmc_toolkit.scheme_d` and serves as an oracle for the compiler |
| scheme D | the package of R1, its interpreter and the independent enumeration oracle (`index_oracle`) |
| DiagramIR | serializable description of a diagram family ([diagram-compiler.md](diagram-compiler.md)) |
| cheap evaluator | evaluator derived from the same IR under a policy that drops the expensive factors |
| NativeKernelIR | typed control-flow form of the delayed-acceptance kernel, from which Fortran is emitted |
| stage 1–6 | the validation steps of the compiler, `scripts/dc_stage1_exact.py` … `scripts/dc_stage6_fepdmc_toolchain.py` |
| R1 grouped toy | finite model that tests whether a grouped measure built from summed window operators can be sampled exactly ([r1-grouped-toy.md](r1-grouped-toy.md)) |

## Records

The files under `benchmarks/`, `validation/` and `provenance/` are kept as
they were produced. Paths and identifiers inside them (for example
`research/…`, `closure/…`, `future-b-…`) refer to the environment of the
original runs, not to this repository.
