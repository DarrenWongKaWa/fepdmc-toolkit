# Tests

```bash
pytest -q
```

No FEP-DMC checkout, material data or Docker needed. The tests in
`test_fortran_backend.py` and `test_fortran_primitives.py` need `gfortran`
(override with `FEPDMC_TOOLKIT_GFORTRAN`) and are skipped without it. Six
more tests are skipped unless the native pieces are present:

| tests | need |
|---|---|
| `test_prepare_upstream.py` (4) | `FEPDMC_TOOLKIT_FEP_DMC` pointing at a clone of the pinned FEP-DMC commit |
| `test_dc_fep_dmc_validate.py` (2) | the same clone, and the Docker image named by `FEPDMC_TOOLKIT_DC_DOCKER` (default `fepdmc-toolkit-u22:local`) |

`tests/_synthetic_runs.py` writes deterministic stand-ins for run directories;
`tests/data/analysis_reference.json` holds the outputs of the analyses on them.
What each group checks is listed in
[docs/validation.md](../docs/validation.md).
