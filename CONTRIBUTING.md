# Contributing

## Set up and check

```bash
python -m venv .venv
.venv/bin/pip install -e ".[dev]" ruff==0.16.9
.venv/bin/ruff check .
.venv/bin/ruff format --check .
.venv/bin/pytest -q
```

`gfortran` is needed for the tests that compile generated Fortran; without
it they are skipped. The continuous integration runs the same commands on
Python 3.11 to 3.14, the tests from an unpacked source distribution, the
validation stages of the diagram compiler and the patches on a clone of the
pinned FEP-DMC commit.

To run the tests that need the pinned FEP-DMC clone and the Docker image:

```bash
git clone https://github.com/yaoluo/FEP-DMC /tmp/fep-dmc
git -C /tmp/fep-dmc checkout 05d08449cffdbd0dfbbbf5009add5cc887bc754b
docker build -f integration/diagram_compiler/Dockerfile.u22 -t fepdmc-toolkit-u22:local .
FEPDMC_TOOLKIT_FEP_DMC=/tmp/fep-dmc .venv/bin/pytest -q
```

## Files that must not be edited

`validation/SHA256SUMS` lists the validated artifacts: the Fortran that is
patched into FEP-DMC or compiled, the recorded benchmark results, the
upstream pins and the validation records. `tests/test_integrity.py` fails if
one of them changes.

- A change to a listed Fortran file or to the text a patch inserts changes
  the patched FEP-DMC source. It needs a new native validation
  ([docs/validation.md](docs/validation.md)); record the new checksum in the
  same commit.
- `provenance/PREPARED_TREE.json` holds the hashes of the tree that
  `fepdmc-toolkit prepare` produces from the pinned commit.
  `tests/test_prepare_upstream.py` compares against it. A refactoring of the
  patch code must leave these hashes unchanged.
- The files under `benchmarks/` and `validation/` are records of runs. They
  are never edited; a new run adds a new record.

## Conventions

- New code follows the `ruff` configuration in `pyproject.toml`.
- A patch names the exact text it replaces and raises `PatchError` if the
  text is not found exactly once. Patches accept only the pinned commit.
- Command-line tools return 0 on success, 1 if a check fails or a patch is
  refused, and 2 if the input is unusable. Errors go to standard error;
  `compare` and `analyze` print their result as JSON on standard output.
- No statement about speed without a recorded benchmark and its statistics
  ([docs/statistics.md](docs/statistics.md)).
- The run-time switches of the patched binary keep the prefix `FUTUREB_`
  ([docs/glossary.md](docs/glossary.md)).

## Release

The version is defined once, in `src/fepdmc_toolkit/__init__.py`. For a
release, set it there, in `CITATION.cff` and in a new section of
`CHANGELOG.md`, then tag `v<version>`.
