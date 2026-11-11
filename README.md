# fepdmc-toolkit

[![ci](https://github.com/DarrenWongKaWa/fepdmc-toolkit/actions/workflows/ci.yml/badge.svg)](https://github.com/DarrenWongKaWa/fepdmc-toolkit/actions/workflows/ci.yml)
[![native-validation](https://github.com/DarrenWongKaWa/fepdmc-toolkit/actions/workflows/native-validation.yml/badge.svg)](https://github.com/DarrenWongKaWa/fepdmc-toolkit/actions/workflows/native-validation.yml)

Source patches, run tools, statistics and a diagram compiler for
[FEP-DMC](https://github.com/yaoluo/FEP-DMC), the first-principles
electron–phonon diagrammatic Monte Carlo code of Luo, Park and Bernardi
([Nat. Phys. 21, 1275 (2025)](https://doi.org/10.1038/s41567-025-02954-1)),
built on Perturbo and Quantum ESPRESSO.

This is an independent toolkit. It is not part of FEP-DMC or Perturbo and
does not redistribute them: every patch is applied to a user-provided
checkout of the public FEP-DMC repository at the pinned commit
`05d08449cffdbd0dfbbbf5009add5cc887bc754b`.

## Components

| Component | What it does | Entry point | Docs |
|---|---|---|---|
| FEP-DMC patches and runner | Opt-in fixes for three reported defects in the native external-phonon updates, exact Monte Carlo moves, diagnostics, a Docker run driver, an end-to-end exactness check, a pooled ratio estimator and analyses of run directories | `fepdmc-toolkit` command, `fepdmc_toolkit.fepdmc` | [docs/fepdmc.md](docs/fepdmc.md) |
| C0 adapter | Source transform that refreshes the phonon frequency after external-phonon sampling in `add_external_ph` | `integration/c0/` | [docs/c0-adapter.md](docs/c0-adapter.md) |
| P1 adapter | Exact-corrected two-stage delayed acceptance for the swap update | `integration/p1/`, `fepdmc_toolkit.p1_da` | [docs/p1.md](docs/p1.md) |
| Statistics | Formation energy from the native estimator; jackknife errors, the HAC stability rule and the bootstrap classification rule for correlated chains | `fepdmc_toolkit.formation_energy`, `fepdmc_toolkit.stats` | [docs/statistics.md](docs/statistics.md) |
| Diagram compiler (experimental) | Fixed-order diagram families → shared DAG → cheap evaluator → delayed-acceptance kernel → NativeKernelIR → generated Fortran | `fepdmc_toolkit.diagram_compiler` | [docs/diagram-compiler.md](docs/diagram-compiler.md) |
| Upstream pin | Fail-closed verification of a FEP-DMC checkout against the pinned commit | `tools/verify_fep_dmc_upstream.py` | [docs/upstream.md](docs/upstream.md) |

The short names used in the code and the records (C0, P1, R1, `FUTUREB_…`,
B-best, …) are explained in [docs/glossary.md](docs/glossary.md).

## Install

Python 3.11 or newer.

```bash
git clone https://github.com/DarrenWongKaWa/fepdmc-toolkit
cd fepdmc-toolkit
python -m venv .venv
.venv/bin/pip install -e ".[dev]"
.venv/bin/pytest -q
```

The tests need no FEP-DMC checkout, no material data and no Docker. The
tests that compile generated Fortran need `gfortran` (set
`FEPDMC_TOOLKIT_GFORTRAN` to choose another compiler) and are skipped
without it. Six tests that need a clone of the pinned FEP-DMC commit or the
Docker image are skipped otherwise ([tests/README.md](tests/README.md)).

## Quick start

Patch a FEP-DMC checkout, build it, run chains and check exactness (details
and all switches in [docs/fepdmc.md](docs/fepdmc.md)):

```bash
git clone https://github.com/yaoluo/FEP-DMC
git -C FEP-DMC checkout 05d08449cffdbd0dfbbbf5009add5cc887bc754b
python tools/verify_fep_dmc_upstream.py FEP-DMC

fepdmc-toolkit prepare FEP-DMC/perturbo-fep-dmc --profile fixes --dry-run   # try the patches, write nothing
fepdmc-toolkit prepare FEP-DMC/perturbo-fep-dmc --profile fixes             # or: --profile research
# build pert-src inside your QE 6.5 tree, then:
export FEPDMC_TOOLKIT_WS=<workspace> FEPDMC_TOOLKIT_DATA=<FEP-DMC dataset>
fepdmc-toolkit run prepare-tables lif_hole --alias FEP-DMC/example/pert-Htable/alias.py
fepdmc-toolkit run ez lif_hole 8 --build <build> --wqfix --extfix --extrmfix   # add --dry-run to see the commands
fepdmc-toolkit validate --build <build>
```

`prepare` writes the tree only if every patch applies. A tree in which the
patched text differs, such as the current upstream tip or a tree that is
already patched, is refused and left as it is.

Apply the delayed-acceptance adapter instead ([docs/p1.md](docs/p1.md)):

```bash
python integration/p1/apply_p1.py --tree FEP-DMC --dry-run
python integration/p1/apply_p1.py --tree FEP-DMC
python integration/p1/verify_p1.py --tree FEP-DMC
```

Compile and evaluate a diagram family in Python
([docs/diagram-compiler.md](docs/diagram-compiler.md)):

```bash
python examples/diagram_compiler/quickstart.py
```

## What has been validated

| Item | Level | Record |
|---|---|---|
| External-phonon fixes and exact moves | Native LiF and SrTiO₃ runs; independent kernels agree at maxOrder 7 (`fepdmc-toolkit validate`) | [docs/fepdmc.md](docs/fepdmc.md), `validation/fepdmc/` |
| `fepdmc-toolkit prepare` | The patched tree is reproduced file by file on the pinned commit (both profiles) | `provenance/PREPARED_TREE.json` |
| C0 adapter | Source transform against the pinned commit (exact pre- and post-image hashes) | `provenance/C0_PUBLIC_PATCH.json` |
| P1 adapter | Source transform plus a native build and a short public LiF fixture | [docs/p1.md](docs/p1.md), `validation/p1/` |
| P1 in a LiF benchmark (benchmark build) | Event-level checks of the exact target; wall-time ratio to the native build 0.956, bootstrap interval [0.940, 0.972]; a ≥ 5% gain is not established | [docs/results.md](docs/results.md) |
| Diagram compiler | 224 differential cases against two independent oracles (max. abs. error 2.8 × 10⁻¹⁷); the generated Fortran, compiled with gfortran 11.4 on Ubuntu 22.04, matches the IR interpreter | [docs/diagram-compiler.md](docs/diagram-compiler.md), `validation/diagram_compiler/` |

No component has an established speed-up. The full list of what is and is
not established is in [docs/validation.md](docs/validation.md) and
[docs/limitations.md](docs/limitations.md).

## Repository layout

```
src/fepdmc_toolkit/   Python package (patches, runner, statistics, compiler)
integration/          C0, P1 and diagram-compiler adapters for the pinned FEP-DMC tree
tools/                upstream pin verifier
scripts/              diagram-compiler validation stages 1-6
provenance/           pinned upstream commit, patch identities, prepared-tree hashes
benchmarks/           recorded LiF benchmark results
validation/           validation records and checksums of validated artifacts
data/, examples/      small public fixture and runnable examples
docs/                 documentation (index: docs/README.md)
tests/                test suite
```

## Development

`ruff check .`, `ruff format --check .` and `pytest -q` must pass; the
workflows in `.github/workflows/` run them, the validation stages and the
adapters on the pinned upstream tree. The validated Fortran and the recorded
results must not be edited ([CONTRIBUTING.md](CONTRIBUTING.md)).

## License and attribution

GPL-3.0-or-later. The Fortran sources under `src/fepdmc_toolkit/fortran/`,
`src/fepdmc_toolkit/fepdmc/data/` and `integration/p1/` are derived from or
written against Perturbo/FEP-DMC (GPL-3). FEP-DMC, Perturbo, the
electron–phonon compression and the material datasets belong to their
authors; cite Luo, Park and Bernardi (2025) and Perturbo when you use this
toolkit. See [NOTICE](NOTICE), [docs/upstream.md](docs/upstream.md) and
[CITATION.cff](CITATION.cff).
