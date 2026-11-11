# Recorded benchmark results

Summaries of the LiF benchmark runs, described in
[docs/results.md](../docs/results.md). No material data.

| directory | contents |
|---|---|
| `c5_dev/` | formation energy Q for three builds, six chains each, with per-chain HAC diagnostics |
| `p1_vs_bbest/` | P1 against the native build: Q, wall time, bootstrap classification, counters, binary hashes |
| `c0_c1/` | C0 frequency-refresh acceptance check and the swap-update profile |
| `r1_fixed_order/` | R1 evaluator against optimized CSE |

The files are checked by the tests and listed in `validation/SHA256SUMS`;
do not edit them. They are kept as produced: paths and identifiers inside
them refer to the environment of the original runs
([docs/glossary.md](../docs/glossary.md)).
