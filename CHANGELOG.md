# Changelog

## 1.0.0 — 2026-11-11

First release as `fepdmc-toolkit`, reorganized from an earlier research
codebase. The patched and generated Fortran, the recorded results and the
upstream pins are byte-identical to the validated versions
(`validation/SHA256SUMS`).

Package and command line

- One Python package, `fepdmc_toolkit`, with the FEP-DMC patches and runner,
  the P1 algebra, statistics, the diagram compiler and its reference oracles.
- Command `fepdmc-toolkit` with `prepare`, `run`, `validate`, `compare`,
  `analyze` and `--version`.
- `prepare` tries the patches on a copy and writes the tree only if all of
  them apply; `--dry-run` lists the files. A refused patch raises
  `PatchError` and leaves the tree unchanged. The result on the pinned
  commit is recorded in `provenance/PREPARED_TREE.json`.
- `run` checks for Docker, the dataset, the binary and the tables before it
  starts a container; `--dry-run` prints the docker commands. The site is
  configured with `FEPDMC_TOOLKIT_WS`, `FEPDMC_TOOLKIT_DATA` and
  `FEPDMC_TOOLKIT_IMAGE`.
- `analyze <analysis>` runs the analyses of run directories with a common
  interface (`--list`, `--help`); each analysis also returns its result as a
  dict. `summarize_sign` reads the table written by `collect_sign`.
- `compare` and `validate` report missing or unfinished chains instead of
  failing with an arithmetic error.

Adapters and validation

- C0 and P1 source adapters and the upstream verifier for FEP-DMC commit
  `05d08449`.
- Diagram-compiler validation scripts `scripts/dc_stage1_exact.py` …
  `dc_stage6_fepdmc_toolchain.py`. Each stage writes `summary.json` and exits
  with status 1 if a check fails.
- Tests of the analyses on deterministic stand-ins for run directories, and
  of `prepare` on a clone of the pinned commit.
- Continuous integration on Python 3.11 to 3.14 with `ruff`, the tests, the
  source distribution, the validation stages and the pinned upstream tree.

Documentation

- Documentation rewritten, with a glossary of the short names
  (`docs/glossary.md`). The upstream status of the reported FEP-DMC defects
  is recorded: the frequency fix was merged upstream in `8ad55d1`; the two
  `remove_external_ph` items are disputed.
