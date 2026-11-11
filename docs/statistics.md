# Statistics

The protocol used for the LiF benchmark in [results.md](results.md), and the
helpers in `fepdmc_toolkit.formation_energy` and `fepdmc_toolkit.stats`:
jackknife standard errors, the HAC stability rule and the bootstrap
classification rule. The package does not include a HAC estimator or a
bootstrap engine; the rules take their per-lag standard errors and interval
end points as inputs.

## Formation energy

The native EZ estimator gives a ratio of block sums. The polaron formation
energy is that ratio minus the bare-band energy:

```
Q = Re( sum_b N_b / tauMax_b  /  sum_b D_b ) - E_bare
```

For the LiF-electron setting used here (Method 0, 20³ phonon grid, rank 20),
`E_bare = 9.35487318746596053 eV` (`E_BARE_LIF_METHOD0`). The unshifted ratio,
about +9.105 eV, is not Q. The shift does not change method-to-method
variances or timing ratios, but it changes the meaning of the number and the
scale of a relative error (1% of |Q| is about 2.5 meV).

## Per-chain HAC

Bartlett HAC standard errors on production blocks at lags 2, 4 and 8. A chain
fails the stability rule if max SE / min SE > 1.25 across the three lags
(`hac_pass`). Failing chains are reported, not dropped, and the lags are fixed
before looking at the data.

## Across-chain jackknife

Independent replicate chains, delete-one jackknife SE of Q
(`mean_jackknife_se`, `ratio_of_sums_jackknife_se`). A point 95% half-width
uses Student t with n − 1 degrees of freedom (2.57058 for six chains).
Relative figures use |Q| as the denominator. These are statistical
diagnostics, not a total uncertainty including projection time, order cut-off
and numerical remainder.

## Timing and cost

`T` is the host wall time of the same Docker run that produced the blocks.
The cost at equal precision is diagnosed by `K = T · Var_JK(Q)`.

## Bootstrap classification of a speed-up

2000 replicate-index resamples, joint across the two arms. With a threshold
of 0.95 on the wall-time ratio (a ≥ 5% improvement), `classify_p1_wall_ratio`
returns:

| rule | classification |
|---|---|
| bootstrap 97.5% < 0.95 and HAC passes in both arms | `P1_GAIN_ESTABLISHED` |
| bootstrap 2.5% > 0.95 | `P1_NO_PRACTICAL_GAIN` |
| the interval crosses 0.95, or HAC fails | `P1_UNRESOLVED_WITHIN_BUDGET` |

An interval that contains 1 is not a proof of equality, and an interval that
crosses the threshold is not a reason to keep sampling until it does not:
fix the budget before the runs.

## Low-sign runs

For FEP-DMC runs with a small average sign, pool the ratio Σn/Σd over chains
instead of averaging per-chain ratios. `fepdmc-toolkit compare` does this and
reports the chain-level variance, an F interval for the variance ratio (needs
scipy), the CPU ratio and the efficiency ([fepdmc.md](fepdmc.md)).
