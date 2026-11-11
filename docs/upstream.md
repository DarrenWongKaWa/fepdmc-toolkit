# Upstream pin, attribution and licensing

## Pinned FEP-DMC commit

| field | value |
|---|---|
| repository | https://github.com/yaoluo/FEP-DMC |
| reference commit | `05d08449cffdbd0dfbbbf5009add5cc887bc754b` (2025-09-24) |
| record | [provenance/UPSTREAM_FEP_DMC.json](../provenance/UPSTREAM_FEP_DMC.json) |
| patch records | [C0](../provenance/C0_PUBLIC_PATCH.json), [P1](../provenance/P1_PUBLIC_PATCH.json), [`prepare`](../provenance/PREPARED_TREE.json): hashes of the files each transform produces from this commit |

The pinned files are `diagMC_JJ_updates.f90` (swap, `add_external_ph`),
`diagMC.f90` (`Gel`, `Dph`), `pert_param.f90` (`zeroTMC`, `DMC_Method`) and
`pert-src/makefile`. Their native Fortran is unchanged from the first public
source upload `0f9b7aa` (2024-09-02) to the pinned commit.

Later upstream commits: `8ad55d1` (2026-10-02) fixes the stale phonon
frequency in `add_external_ph`, the defect that C0 and `--wqfix` address
([yaoluo/FEP-DMC#1](https://github.com/yaoluo/FEP-DMC/issues/1)); `5ec3bd1`
edits the README. The patches in this repository target the pinned commit.

### Verify a checkout

```bash
git clone https://github.com/yaoluo/FEP-DMC
git -C FEP-DMC checkout 05d08449cffdbd0dfbbbf5009add5cc887bc754b
python tools/verify_fep_dmc_upstream.py ./FEP-DMC
```

The tool does not fetch, check out or edit the target. Exit 0 requires the
directory to be its own git top level, HEAD to equal the pinned commit, the
pinned files' SHA256 to match the on-disk bytes, and those files to be clean
in `git status`. A snapshot without its own `.git` may match the hashes and
still report `GIT_UNAVAILABLE` (exit 3): content identity, not full upstream
identity. Hashes are of the bytes GitHub serves (LF line endings);
`core.autocrlf` conversion fails the check. A patched tree fails the check
by design.

### The benchmarked builds

The LiF benchmark builds in [results.md](results.md) were made from a
private working tree of the same upstream code, not from the pinned public
commit; byte equivalence between the two is not established. The binaries
are identified by SHA256 in `benchmarks/`, which are not source revisions.
The Fortran of the timed P1 build is kept under `src/fepdmc_toolkit/fortran/`
(`linear_da_mod.f90` SHA256 of the timed file:
`004b23b66e00f4ca20e2cc683c8e99124f12713549fcc02363e0e0afc37854df`; the
shipped copy adds a fixture-only gate that is inactive in production).

## Attribution

Not part of this toolkit, and not claimed by it:

- the FEP-DMC method, implementation and LiF results: Y. Luo, J. Park and
  M. Bernardi, Nat. Phys. **21**, 1275 (2025),
  doi:[10.1038/s41567-025-02954-1](https://doi.org/10.1038/s41567-025-02954-1);
- Perturbo (GPL-3), https://perturbo-code.github.io/;
- the Luo–Bernardi electron–phonon compression and the matrix-product band
  summation used by the native code;
- the upstream material datasets (epwan files, H tables, SVD data);
- Quantum ESPRESSO 6.5 (GPL), linked by the native build.

| path | origin | license |
|---|---|---|
| `src/fepdmc_toolkit/fortran/linear_da_mod.f90` | FEP-DMC `pert-src` plus the delayed-acceptance helpers of the benchmarked build | GPL-3 |
| `src/fepdmc_toolkit/fortran/update_swap_p1_excerpt.f90` | excerpt of `diagMC_JJ_updates.f90` `update_swap` with the P1 insertion | GPL-3 |
| `src/fepdmc_toolkit/fortran/closure_observer.f90` | timers of the benchmarked build | GPL-3 |
| `src/fepdmc_toolkit/fepdmc/data/bchain.f90`, `rb_window.f90` | this toolkit, written against FEP-DMC `pert-src` interfaces | GPL-3 |
| `src/fepdmc_toolkit/fepdmc/data/mkl_vsl.f90` | this toolkit: xorshift64 shim for the MKL VSL symbols used by FEP-DMC when MKL is absent | GPL-3 |
| `src/fepdmc_toolkit/fepdmc/patches/*.py`, `integration/*` | this toolkit: anchored source transforms; the patched files remain Perturbo-derived | GPL-3 |
| other Python code | this toolkit | GPL-3 |

## Data not included

| item | reason |
|---|---|
| `lif-sp3_epwan.h5` and other epwan HDF5 | upstream material data, hundreds of MiB to GiB, no redistribution grant |
| the Nature Physics PDF | publisher copyright |
| Docker image layers | build environment, not source |
| QE `libpw.a`, `.o`, `.mod`, `perturbo.x` | compiled artifacts; rebuild from upstream |
| Intel compiler runtimes | proprietary |

Obtain FEP-DMC from https://github.com/yaoluo/FEP-DMC and the datasets from
the FEP-DMC authors. When you fork this repository, keep GPL-3 on the
Perturbo-derived Fortran, attribute Luo, Park, Bernardi and Perturbo, and do
not add the datasets or the paper.
