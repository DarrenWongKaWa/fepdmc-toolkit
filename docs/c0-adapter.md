# C0 adapter: phonon-frequency refresh in `add_external_ph`

C0 is a source transform of the pinned FEP-DMC tree. It inserts one line
and writes one file:

`perturbo-fep-dmc/pert-src/diagMC_JJ_updates.f90`

Upstream fixed the same defect in commit `8ad55d1` (2026-10-02,
[yaoluo/FEP-DMC#1](https://github.com/yaoluo/FEP-DMC/issues/1)). C0 is
needed only on the pinned commit `05d08449`, which the other patches in this
repository target.

## What C0 changes

In `add_external_ph`, `sample_q_omp_int` writes a new `vn1%i_q` and does not
set `vn1%wq`. Later the same update copies `vn2%wq = vn1%wq`, samples
imaginary time with `vn1%wq(vn1%nu)`, and evaluates `Dph` with those
frequencies. The sibling update `add_ph` already calls
`cal_wq_int( vn1%i_q, vn1%wq )` after its sampler. C0 inserts that same
native call after the external sampler.

`sample_q_omp_int` (in `diagMC.f90`) samples a grid index into `ixqt`
(`vn1%i_q`) and a proposal probability. It does not compute phonon
frequencies. `cal_wq_int` is the native grid-index frequency routine
(`solve_phonon_fast`, then Ry→eV). Without it, `add_external_ph` keeps
whatever `wq` the recycled vertex slot last held.

## Supported upstream

| field | value |
|---|---|
| repository | https://github.com/yaoluo/FEP-DMC |
| commit | `05d08449cffdbd0dfbbbf5009add5cc887bc754b` |
| preimage SHA256 | `6c97c694f76df18e9adee7b79666b480f1de5b00f30e8bb5ba08510289e43571` |
| postimage SHA256 | `be3a2215ab257e12aaaba1a0a3cca8714a4174c0fee890bb881cd75a32b9a235` |

Machine record: [provenance/C0_PUBLIC_PATCH.json](../provenance/C0_PUBLIC_PATCH.json).

## Use

```bash
git clone https://github.com/yaoluo/FEP-DMC
git -C FEP-DMC checkout 05d08449cffdbd0dfbbbf5009add5cc887bc754b
python tools/verify_fep_dmc_upstream.py ./FEP-DMC     # pristine identity

python integration/c0/apply_c0.py --tree ./FEP-DMC --dry-run   # unified diff, zero writes
python integration/c0/apply_c0.py --tree ./FEP-DMC
python integration/c0/verify_c0.py --tree ./FEP-DMC            # C0_STATE=PRISTINE|APPLIED|UNKNOWN
```

The tools do not fetch, check out or reset your tree. Identity requires the
directory to be its own git top level at the pinned commit.

| state | apply |
|---|---|
| PRISTINE | insert exactly one `cal_wq_int` after the external sampler |
| APPLIED | print `ALREADY_APPLIED`, exit 2, no write |
| UNKNOWN | refuse, exit 1, no write |

UNKNOWN includes a wrong commit, dirty unrelated pinned files, a missing or
duplicated anchor, partial edits, and snapshots without git identity. A
C0-patched tree does not pass `tools/verify_fep_dmc_upstream.py`, which
stays a pristine-identity check. C0 composes with the P1 adapter in either
order ([p1.md](p1.md)).

Revert from the FEP-DMC checkout (the adapter never runs git itself):

```bash
git restore perturbo-fep-dmc/pert-src/diagMC_JJ_updates.f90
```

## Scope

C0 is validated as a source transform (exact pre- and post-image hashes).
The adapter itself does not compile `perturbo.x` or run a material
calculation. The same one-line fix, as the `--wqfix` switch of
`fepdmc-toolkit`, was used in the native runs described in
[fepdmc.md](fepdmc.md).
