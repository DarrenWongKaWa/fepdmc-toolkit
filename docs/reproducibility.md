# Reproducibility

## Without the native stack

From a clone, with Python 3.11 or newer:

```bash
pip install -e ".[dev]"
pytest -q                                          # full test suite
python examples/lif_small_fixture/recompute_q.py   # recompute the recorded LiF summaries
python examples/diagram_compiler/quickstart.py     # compile and evaluate a diagram family
```

The diagram-compiler validation stages 1–5 run without FEP-DMC and write
their records to `build/validation/` (stage 5 needs gfortran). Each stage
prints PASS or FAIL, writes `summary.json` with its checks and exits with
status 1 on a failure:

```bash
python scripts/dc_stage1_exact.py            # 224 differential cases against R1 and the index oracle
python scripts/dc_stage2_cheap.py
python scripts/dc_stage3_da_kernel.py
python scripts/dc_stage4_native_ir.py
python scripts/dc_stage5_fortran_backend.py
```

Stage 6 needs Docker and a clone of the pinned FEP-DMC commit; see
[integration/diagram_compiler/README.md](../integration/diagram_compiler/README.md).

The source adapters can be checked against the real upstream tree without
compiling:

```bash
git clone https://github.com/yaoluo/FEP-DMC /tmp/fep
git -C /tmp/fep checkout 05d08449cffdbd0dfbbbf5009add5cc887bc754b
python tools/verify_fep_dmc_upstream.py /tmp/fep
python integration/c0/apply_c0.py --tree /tmp/fep && python integration/c0/verify_c0.py --tree /tmp/fep
git -C /tmp/fep restore perturbo-fep-dmc/pert-src/diagMC_JJ_updates.f90
python integration/p1/apply_p1.py --tree /tmp/fep && python integration/p1/verify_p1.py --tree /tmp/fep
git -C /tmp/fep checkout -- . && git -C /tmp/fep clean -fdq
FEPDMC_TOOLKIT_FEP_DMC=/tmp/fep pytest -q tests/test_prepare_upstream.py   # prepare, both profiles
```

The CI workflows in `.github/workflows/` run these steps.

## Native runs

| piece | source | included |
|---|---|---|
| Quantum ESPRESSO 6.5 | https://gitlab.com/QEF/q-e/-/tags/qe-6.5 | no |
| FEP-DMC / Perturbo | https://github.com/yaoluo/FEP-DMC at the pinned commit | patches only |
| LiF `lif-sp3_epwan.h5` and other datasets | FEP-DMC authors | no |
| `gkq-20` tables | built from the epwan files (`fepdmc-toolkit run prepare-tables`) | no |
| Fortran compiler with HDF5 (Fortran 2003) | vendor or distribution | no |
| Docker image with the QE 6.5 toolchain (`r5p0-env:ubuntu2004` in the recorded runs) | local build | no |

The recorded LiF benchmarks used the Docker image digest
`sha256:2082c5082f050fafc807ad1488e3498fe51e0c7951e4ddcf818b3ed39de5f7e1`
(linux/amd64), Intel Fortran, `OMP_NUM_THREADS=1` and one chain per 6 GiB
container slot. Seeds: master 9143600 (spawn 0 = B-best, 1 = P1, six
children each), bootstrap 9143610, `Nmcmc(1e4) = 200`.

The native validation of the P1 adapter used Ubuntu 22.04 with gfortran 11.4,
OpenMPI and the distribution HDF5 instead
([validation/p1/NATIVE_VALIDATION.json](../validation/p1/NATIVE_VALIDATION.json)).
