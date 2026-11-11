# FEP-DMC patches and runner (`fepdmc_toolkit.fepdmc`)

Source patches, a Docker run driver and estimators on top of upstream
FEP-DMC (Luo, Park, Bernardi), pinned at
`05d08449cffdbd0dfbbbf5009add5cc887bc754b`
([provenance/UPSTREAM_FEP_DMC.json](../provenance/UPSTREAM_FEP_DMC.json)).
Three parts:

- fixes for three defects in the native external-phonon updates, all opt-in at run time;
- exact Monte Carlo moves that the native code does not have;
- diagnostics and estimators for the low-sign regime.

The patches are source transforms of a user-provided tree plus Python
drivers. Nothing here is a compiled plugin, and no move has an established
efficiency gain.

## Install and quick start

```bash
pip install .                   # from a checkout; installs the fepdmc-toolkit command

# 1. patch a pristine checkout at the pin
#    profile "research" = fixes + moves + diagnostics; "fixes" = fixes only
fepdmc-toolkit prepare <path>/perturbo-fep-dmc --profile research --dry-run   # optional trial
fepdmc-toolkit prepare <path>/perturbo-fep-dmc --profile research

# 2. build pert-src inside the QE 6.5 tree (gfortran, serial make)
cd <ws>/build/q-e-qe-6.5/<build>/pert-src && make

# 3. run chains in Docker
export FEPDMC_TOOLKIT_WS=<ws>  FEPDMC_TOOLKIT_DATA=<FEP-DMC dataset dir>
fepdmc-toolkit run prepare-tables lif_hole --alias <FEP-DMC>/example/pert-Htable/alias.py
fepdmc-toolkit run ez lif_hole 8 --build <build> --temp 500 --maxorder 7 \
    --wqfix --extfix --extrmfix --tag my_runs

# 4. end-to-end exactness check (independent kernels must agree)
fepdmc-toolkit validate --build <build>

# 5. compare two sets of runs with the pooled estimator
fepdmc-toolkit compare <runs_A> <runs_C> <material> <E_bare> <tau_max>

# 6. further analyses of run directories
fepdmc-toolkit analyze --list
```

Each patch names the exact text it edits and refuses a tree in which that
text does not occur exactly once. `prepare` first applies the patches to a
temporary copy and writes the tree only if all of them apply, so a refused
tree (another commit, or a tree that is already patched) is left unchanged;
`--dry-run` stops after the trial and lists the files. Applied to a clean
pin, `prepare` reproduces the sources of the validated research build byte
for byte. The hashes of the files it writes, for both profiles, are in
[provenance/PREPARED_TREE.json](../provenance/PREPARED_TREE.json) and are
checked by `tests/test_prepare_upstream.py`.

`run` checks that Docker, the dataset, the binary and the tables exist
before it starts a container. `run … --dry-run` prints the docker commands
and writes nothing.

### Configuration

| flag | environment variable | default |
|---|---|---|
| `--ws` | `FEPDMC_TOOLKIT_WS` | `./fepdmc_ws`, mounted at `/work` |
| `--data` | `FEPDMC_TOOLKIT_DATA` | required; mounted read-only at `/data` |
| `--image` | `FEPDMC_TOOLKIT_IMAGE` | `r5p0-env:ubuntu2004` (QE 6.5 toolchain) |
| `--build` | — | `perturbo-fep-dmc-pkg`, i.e. `<ws>/build/q-e-qe-6.5/<build>/pert-src/perturbo.x` |
| `--tables` | — | `<ws>/tables/<material>/gkq-20` |

Runs land in `<ws>/<tag>/<material>/chainNN` with `run_meta.json` (seed,
switches, binary). Materials: `lif_elec`, `lif_hole`, `sto`, `anatase`,
`anatase_b3`. `tabulate-H` needs a build made without
`-fmax-stack-var-size=1`; pass it with `--tables-build`.

The switches of the patched binary are environment variables with the
prefix `FUTUREB_` ("Future B" was the working name of the project). The
prefix is part of the patched Fortran, which is kept byte-identical to the
validated builds. `run` sets the switches inside the container from its
flags; set them yourself only if you start `perturbo.x` by hand.

## Fixes for the native external-phonon updates

| switch (`run` flag) | patch | defect | upstream status |
|---|---|---|---|
| `FUTUREB_WQFIX=1` (`--wqfix`) | [wqfix](../src/fepdmc_toolkit/fepdmc/patches/wqfix.py) | `add_external_ph` never calls `cal_wq_int`, so external pairs keep the stale phonon frequency of the recycled vertex slot. The C0 adapter ([c0-adapter.md](c0-adapter.md)) fixes the same line unconditionally | fixed upstream in `8ad55d1` (2026-10-02) |
| `FUTUREB_EXTRMFIX=1` (`--extrmfix`) | [extrmfix](../src/fepdmc_toolkit/fepdmc/patches/extrmfix.py) | `remove_external_ph` builds the pair-less reference trace at momentum k − q instead of k; multiband only | reported; disputed by the maintainer, discussion open |
| `FUTUREB_EXTFIX=1` (`--extfix`) | [extfix](../src/fepdmc_toolkit/fepdmc/patches/extfix.py) | `remove_external_ph` evaluates the reverse-add density without a range check, so pairs outside add's τ support are removable | reported; disputed by the maintainer, discussion open |

All three were reported in
[yaoluo/FEP-DMC#1](https://github.com/yaoluo/FEP-DMC/issues/1). The patches
target the pinned commit `05d08449`, which predates the upstream frequency
fix; check out that commit to use them.

`--wqfix` makes the frequency fix switchable for A/B runs against the
unfixed native code. On a tree with C0 applied, `prepare` still works and
`--wqfix` is redundant. C0 refuses a tree that `prepare` has already
patched.

**Effect of the fixes.** Use all three for the target defined by the
fixed updates.

| case | before the removal fixes | all fixes |
|---|---|---|
| LiF-hole, T = 500 K, maxOrder 7 (total energy E) | −0.87311 ± 0.00074 (wq fix only, 4 chains) | −0.87630 ± 0.00033 (Δ = −0.0032 ± 0.0008, about 4σ; 2.4% of Q) |
| LiF-hole, T = 500 K, full order | ⟨s⟩ 0.090 | ⟨s⟩ 0.111, about 3 fewer external pairs, Q +0.095 eV (1.7σ) |
| LiF-electron and STO at β = 232 (the published cases) | — | unchanged within 1–3 meV |

The published LiF-hole case at β = 232 is three-band and could be affected.
There the average sign cannot be resolved, so the effect is not quantified.
Records: [external_pair_validation.json](../validation/fepdmc/external_pair_validation.json),
[lif_hole500_allfix_vs_wqfix.json](../validation/fepdmc/lif_hole500_allfix_vs_wqfix.json).

## Recommended configuration

| goal | switches |
|---|---|
| native target with all three fixes | `--wqfix --extfix --extrmfix` |
| where external pairs are abundant (low sign, multiband, e.g. LiF-hole at T = 500 K) | the above plus `--extar 0.02` (anything in 0.02–0.1 reaches the same equilibrium) |

**Why the external move is needed there.** From the empty start, native
equilibrates the external-pair count very slowly (τ_int ≳ 1800
measurements). A standard 2×10⁶-step run sits in a quasi-stationary state
with about 9 pairs instead of about 11.5. That biases Q by about 10% and
⟨s⟩ by about 2×. The ordered external move reaches the equilibrium and is
exact (maxOrder 7, 21, 51). p = 0.02, 0.05 and 0.1 all give the same
Q (−1.77 to −1.80) and pair count (11.9–12.0)
([lif_hole500_ext_pscan.json](../validation/fepdmc/lif_hole500_ext_pscan.json)).
CPU time is about 1.9× native at every p. That is mostly physical: the
correct equilibrium has larger diagrams, so native moves cost more. It is
not move overhead, and the move does not reduce per-step variance.

## Exact moves (inert unless enabled)

| switch (`run` flag) | move |
|---|---|
| `FUTUREB_CHQ`, `FUTUREB_CHQ_P` (`--chq P`) | change-q: new momentum for an internal line, span shifted |
| `FUTUREB_CHQ_EXT` (`--chq-ext`) | change-q also for boundary-wrapping external pairs |
| `FUTUREB_ANY`, `FUTUREB_ANY_P` (`--any P`) | general-span add/remove of internal lines (native: span 1 only) |
| `FUTUREB_EXTAR`, `FUTUREB_EXTAR_P` (`--extar P`) | external-pair add/remove at any position, respecting native's ordering (every head-attached external vertex precedes every tail-attached one) |
| `FUTUREB_MV_UNTIL=N` (`--mv-until N`) | apply the moves only for the first N steps (burn-in) |
| `FUTUREB_BCHAIN` (+ `_GROUP`, `_ADDREM`, `_BLOCK`, `_SRULE`, `_REFRESH`, `_PAR`) (`--bchain`, `--block`, `--srule`) | grouped-measure chain B (phonon modes folded over non-crossing lines) |

Moves and chain B need the `research` profile. With the `fixes` profile
only the three fix switches exist.

**Validation.** Every move passes a round trip (C → C′ → C restores log w
and the proposal probabilities to ≤ 4 × 10⁻¹⁵) and per-chain momentum,
gauge and sign self-checks. At maxOrder 7 every exact kernel reproduces the
fixed native code (`fepdmc-toolkit validate`).

**Efficiency.** No move has a demonstrated net gain for Q (net efficiency
relative to the fixed native code: chain B ≈ 0.33, mode Rao–Blackwellization
0.77, change-q 0.71 and 0.86, general add/remove 0.51 at full order and
≈ 1.1 at maxOrder 31, moves in burn-in only 0.70). The moves do mix the
slow variables faster: τ_int of the external-pair count drops from ≳ 1800
to 56, and the chain-to-chain sign spread shrinks 2–6×. For the external
move this corrects the native equilibration bias described above. The run
records are in [validation/fepdmc/](../validation/fepdmc/README.md).

## Diagnostics

- `FUTUREB_SEED` reproducible chains (`run` sets `--seed0 + 7919·i`).
- `FUTUREB_RB`, `FUTUREB_MODES*`, `FUTUREB_QGROUP`, `FUTUREB_MODE_RB`,
  `FUTUREB_MODE_SVD`: measurement-side grouping and mode diagnostics.
- `FUTUREB_CHQ_RT`: round-trip check.
- `FUTUREB_HERM_TEST`: hermiticity of the vertex tables.
- `FUTUREB_EXTAR_OUTER`, `_RPOS`, `_RTAU`, `_MIMIC`, `_UNORDERED`:
  restricted and replica variants of the external move, used to isolate the
  reported defects.
- `FUTUREB_PNU_FRO`: gauge-invariant mode proposal.

Traces: `chq_trace.dat` (numerator, sign, order, external pairs, head |k|²,
internal lines) and `chq_summary.dat` (acceptances and self-checks).

## Estimators

Use `fepdmc-toolkit compare` (Python: `fepdmc_toolkit.fepdmc.pooled_side`). It
pools the ratio Σn/Σd over chains and gives the chain-level variance, an F
interval (needs scipy), the CPU ratio and the efficiency. In low-sign runs,
do **not** average per-chain ratios. Do not trust a variance ratio from
fewer than several batches of about 16 chains: the same method's variance
varied about 5× between batches.

Further analyses of run directories run through `fepdmc-toolkit analyze`
(`--list` names them, `<analysis> --help` shows the arguments). Each prints
JSON and is also a Python function, for example
`fepdmc_toolkit.fepdmc.analysis.analyze_rb.analyze(runs, material)`.

| analysis | reads | result |
|---|---|---|
| `analyze_rb` | `rb_trace.dat` | raw against window Rao–Blackwell estimator; sign of the tiled grouped measure |
| `analyze_modes`, `analyze_q`, `analyze_svd` | `mode_trace.dat`, `q_trace.dat`, `mode_svd_trace.dat` | mode, momentum and SVD-basis grouping diagnostics |
| `analyze_mode_rb` | `mode_rb_trace.dat` | raw against mode Rao–Blackwell estimator on the same chains |
| `estimate_groups` | `lines_dump.dat` | offline estimate of wider mode groups |
| `decompose_native` | `rb_trace.dat` | share of the variance of Q that is not due to the sign |
| `compare_bchain`, `compare_between` | `bchain_trace.dat`, `stdout.log` | native chain A against grouped chain B |
| `summarize_materials` | `stdout.log`, `run_meta.json` | average sign and Q per material |
| `collect_sign`, `summarize_sign` | `sign.dat-*` | table of the average sign of existing runs, and its summary per run group |

`analyze_mode_rb` and `compare_between` need scipy. The `*_trace.dat` and
`lines_dump.dat` files are written by the `research` profile with the
matching `run` flags; `stdout.log` and `sign.dat-*` are native output.

## Layout

| path | content |
|---|---|
| `fepdmc/patches/` | `base` (make.sys, build hooks), `bchain`, `wqfix`, `extfix`, `extrmfix`; registry `ORDER`, `PROFILES`; `PatchError` |
| `fepdmc/data/` | `bchain.f90` (moves, chain B, checks), `rb_window.f90` (measurement-side diagnostics), `mkl_vsl.f90` (MKL VSL random-number shim, used when MKL is absent), `make.sys` |
| `fepdmc/prepare.py`, `runner.py`, `pooled.py`, `validate.py`, `cli.py` | the `fepdmc-toolkit` commands |
| `fepdmc/analysis/` | the analyses behind `fepdmc-toolkit analyze` |
| `tests/test_fepdmc_tools.py`, `tests/test_analysis.py`, `tests/test_prepare_upstream.py` | estimator, registry, CLI and configuration tests; the analyses on stand-ins for run directories; `prepare` on a clone of the pinned commit |

The patched Fortran is GPL-3, derived from Perturbo.
