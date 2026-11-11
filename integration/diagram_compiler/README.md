# Native validation of the diagram compiler (stage 6)

Compiles the generated two-band n = 1 kernel with gfortran in a Docker image
built from `Dockerfile.u22` (Ubuntu 22.04) and checks it against the IR
interpreter. The pinned FEP-DMC tree is a witness only: it is verified to be
pristine and is not modified.

```sh
docker build -f integration/diagram_compiler/Dockerfile.u22 -t fepdmc-toolkit-u22:local .
git clone https://github.com/yaoluo/FEP-DMC /tmp/fep-dmc
git -C /tmp/fep-dmc checkout 05d08449cffdbd0dfbbbf5009add5cc887bc754b
FEPDMC_TOOLKIT_FEP_DMC=/tmp/fep-dmc FEPDMC_TOOLKIT_DC_DOCKER=fepdmc-toolkit-u22:local \
  python scripts/dc_stage6_fepdmc_toolchain.py --fep-dmc /tmp/fep-dmc --image fepdmc-toolkit-u22:local
```

There is no default FEP-DMC path; the script refuses to run without
`--fep-dmc` or `FEPDMC_TOOLKIT_FEP_DMC`, and without the image. It prints
PASS or FAIL and exits with status 1 if a check fails. With the same two
environment variables set, `pytest tests/test_dc_fep_dmc_validate.py` runs
the native tests instead of skipping them.
