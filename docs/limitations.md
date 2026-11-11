# Limitations

## General

- **No established speed-up.** Neither P1, the exact moves nor the diagram
  compiler has a demonstrated net efficiency gain.
- **Pinned upstream.** All patches target FEP-DMC commit `05d08449`. The
  upstream tip has since fixed the stale-frequency defect (`8ad55d1`), so
  the patches and the verifier refuse it.
- **Native stack not included.** Rebuilding `perturbo.x` needs Quantum
  ESPRESSO 6.5, Perturbo/FEP-DMC, a Fortran compiler with HDF5, and the
  material datasets ([reproducibility.md](reproducibility.md)).

## FEP-DMC patches

- Two of the three reported defects (`EXTRMFIX`, `EXTFIX`) are disputed by
  the upstream maintainer; the discussion in
  [yaoluo/FEP-DMC#1](https://github.com/yaoluo/FEP-DMC/issues/1) is open.
- The effect of the fixes on the published three-band LiF-hole case at
  β = 232 is not quantified: the average sign there cannot be resolved.
- Variance ratios between methods varied about 5× between batches of the
  same method; compare methods only over several batches of about 16 chains.

## P1

- Operating domain: JJ `update_swap`, `DMC_Method = 0`, `dmc_band = 1`,
  `sample_gt = .false.`, `zeroTMC = .true.`. Everything else is refused.
- The adapter was validated natively on a short 1000 K fixture only, not on
  production statistics. The LiF benchmark used an earlier build of the same
  algorithm.
- In the benchmark, HAC failed on several chains, so the statistical
  efficiency is unresolved, and a wall-time gain of at least 5% (the
  threshold fixed beforehand) is not established; the point estimate is 4.4%.
- The benchmarked build still contains no-op observer calls; their overhead
  was not subtracted from the wall time. Docker wall time includes container
  start in both arms.
- A full recomputation of Re(M) for the reverse move was not added to the
  production path; the reverse check in the fixture is of the P1 score.

## Statistics and results

- All LiF numbers are for a fixed discrete Hamiltonian (Method 0, 20³, rank
  20, finite projection time and order cap). Projection, order and numerical
  remainders were not closed. They are not thermodynamic-limit, infinite-order
  or 1%-certified ground-state values.

## Diagram compiler

- Two hard-coded families: scalar Holstein (validated n = 1–4) and Hermitian
  two-band (validated n = 1–3). Orders up to 6 are admitted as a resource
  bound but not validated. No variable-order compilation, no topology
  sampler, no general material provider.
- Pairing enumeration is factorial; DAG sharing reduces repeated execution,
  not enumeration.
- Under the chosen cheap policy, the scalar cheap evaluator equals the exact
  one, so it is not cheaper there; the two-band cheap evaluator ignores
  `delta`, `gap`, the σ_z vertices and the trace.
- The compiler's delayed-acceptance kernel targets the positive real diagram
  sum `F`; it is not the P1 target `abs(Re M) · P_kchange`, and production
  FEP-DMC does not call the generated kernel.
- Native validation used one compiler, gfortran 11.4 on Ubuntu 22.04
  (aarch64 and x86_64).
