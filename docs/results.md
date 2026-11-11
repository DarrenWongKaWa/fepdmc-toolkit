# Benchmark results

Recorded results of the LiF benchmarks run with this code. The numbers are
in `benchmarks/` and are checked by the tests; `examples/lif_small_fixture/recompute_q.py`
recomputes the summaries. The runs themselves need the native stack
described in [reproducibility.md](reproducibility.md).

All LiF numbers are for one fixed discrete Hamiltonian: LiF-electron polaron,
Luo EZ estimator, Method 0, 20³ phonon grid, rank 20, finite projection time
and order cap. They are algorithm benchmarks, not converged ground-state
values. Q is the formation energy defined in [statistics.md](statistics.md).

## Formation energy, three builds

Six independent chains per build, N_mc = 200 × 10⁴ attempts each. B0 and
B-best are native builds without delayed acceptance. The P1 row is an earlier
P1 build with its event-level fixture active, so its timing is not used.

| build | Q (eV) | JK SE (meV) | rel. JK SE | 95% half-width (meV) | rel. 95% half-width |
|---|---:|---:|---:|---:|---:|
| B0 | −0.250253260 | 2.368 | 0.946% | 6.088 | 2.433% |
| B-best | −0.247842001 | 3.712 | 1.498% | 9.541 | 3.850% |
| P1 (fixture build) | −0.249583768 | 2.343 | 0.939% | 6.023 | 2.413% |

HAC (Bartlett, lags 2/4/8) failed the stability rule on 5 of 18 chains
(B0_r3, Bbest_r1, Bbest_r2, Bbest_r4, P1_r1). Record:
`benchmarks/c5_dev/`.

Q ≈ −0.25 eV is of the order of Luo et al. Fig. 2(b) at 20³. It is not their
Table 1 value (−0.408 eV) and not the Fan–Migdal/Debye–Waller value
(−0.538 eV).

## P1 against the native build

The P1 build used for timing (binary SHA256 `f9518a21…`) was compiled with
the fixture switch present and run with it off (`P1_FIXTURE=0`,
`LINEAR_DA=on`, `LINEAR_DA_SCORE=prop`). Protocol, threshold and seeds were
fixed before the runs (master seed 9143600, bootstrap seed 9143610).

| build | Q (eV) | JK SE (meV) | rel. SE | wall T (s) | K = T · Var_JK |
|---|---:|---:|---:|---:|---:|
| B-best | −0.251403 | 2.077 | 0.826% | 32.721 | 1.412 × 10⁻⁴ |
| P1 | −0.246614 | 3.350 | 1.358% | 31.293 | 3.511 × 10⁻⁴ |

- Wall-time ratio T_P1 / T_B-best = 0.956, bootstrap 95% interval
  [0.940, 0.972]. The interval crosses the 0.95 threshold, so a ≥ 5% gain is
  not established (`P1_UNRESOLVED_WITHIN_BUDGET`).
- K_P1 / K_B-best = 2.49, bootstrap interval [0.16, 21.5]. The interval is
  too wide to conclude anything, and HAC failed in both arms (B-best chains
  0 and 3, P1 chains 0, 2 and 3). This is not evidence that P1 is slower at
  equal precision.
- Counters over the six P1 chains: 1 395 968 eligible swaps; 463 424 stage-1
  rejects (33.2%); 932 544 reached the exact stage; 30 600 stage-2 rejects;
  901 944 accepts. The identity rejects + accepts = eligible holds.

Record: `benchmarks/p1_vs_bbest/`. The Fortran of the timed build is kept
under `src/fepdmc_toolkit/fortran/`; the public adapter in
`integration/p1/` implements the same algorithm on the pinned FEP-DMC tree
([p1.md](p1.md)).

## Profiling of the native swap update

On LiF, the swap update takes about 21% of the exclusive Monte Carlo time
(`f_swap = 0.212`, 16 profiled runs); environment, propagator and
electron–phonon evaluation remain large costs. Full-detail timers add
about 41% overhead, sparse (1/16) sampling about 3%. Record:
`benchmarks/c0_c1/`.

## R1 fixed-order evaluator against common-subexpression elimination

On a scalar one-dimensional Holstein model at fixed order, compiled once, the
recursive R1 evaluator and an optimized common-subexpression (CSE) baseline
take the same total time: ratio 0.9992 overall, 0.9863–1.0034 across tiers.
Correctness and reuse were validated; no additional speed-up over CSE was
found. Record: `benchmarks/r1_fixed_order/`, metadata in
`fepdmc_toolkit.r1_fixed_order`.

## Not pursued

- A 1% certified ground-state formation energy.
- Learned proposals or neural surrogates.
- A grouped (R1-type) Monte Carlo consumer on a real material; see
  [roadmap.md](roadmap.md).
