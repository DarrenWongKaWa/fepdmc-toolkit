# R1 grouped toy: a finite, enumerable grouped model with two measures

Derivation behind `fepdmc_toolkit.r1_grouped`: can a grouped measure, whose
group weight is a product of summed window operators of the kind R1 computes,
be sampled exactly? No material file is read. Section numbers are cited by
`fepdmc_toolkit.r1_grouped` (`grouped_toy.py`,
`kernels.py` and `run_checks.py`).

## 1. Configuration space and signed weight

Fix a time grid `{1..G}`, momentum torus `Z_L`, maximum order `N`, external
momentum `k_ext` and time step `Δ`. A configuration is

```
C = (τ, P, q):  τ_1 < … < τ_2n ⊂ {1..G},  P a perfect matching of positions 0..2n-1,
                q ∈ Z_L^n indexed by first-opening rank of the lines of P.
```

`Ω` is the finite union over `n = 0..N`. With `τ_0 = 0`, `τ_{2n+1} = G+1`,
segment `j` (ending at vertex `j`) carries `k_j = k_ext − Σ_{(a,b)∈P, a<j≤b} q_(a,b)`.
The weight is

```
D(C) = tr[ U(k_2n, τ_2n+1 − τ_2n) V_2n−1 ⋯ V_0 U(k_0, τ_1) ] · Π_{(a,b)∈P} e^{−ω Δ (τ_b − τ_a)},
U(k, d) = exp(−H(k) Δ d),  H(k) = [[ξ_k, δ], [δ, ξ_k + gap]],  ξ_k = 2t(1 − cos 2πk/L).
```

Vertex matrices depend on the incoming electron momentum (as `g_mnν(k,q)` does):

```
emit  (k → k−q):  M_e(k,q) = (g/√L) (ζ σ_z + η (cos φ σ_x + sin φ σ_y)),  φ = θ_0 + 2πk/L + 0.9 q
absorb (k → k+q): M_a(k,q) = M_e(k+q, q)^†
```

`η = 0, ζ = 1` is exactly the Diagram Compiler two-band family
(`(g/√L) σ_z`, same `H`, Einstein `D`). For `η ≠ 0`, `D(C)` is complex.

The **physical object** is the signed expectation of an observable `O`,

```
⟨O⟩ = Σ_C Re D(C) O(C) / Σ_C Re D(C).
```

Native Luo EZ DiagMC samples `|Re D|` and reweights by the sign. That is
the reference measure **A** below.

## 2. Windows and re-pairing

For `0 ≤ v ≤ 2n−4`, the window `W_v = {v, v+1, v+2, v+3}` is **eligible**
if exactly two lines of `P` touch `W_v` and both lie entirely inside it.
Let `r ∈ {sequential (01)(23), nested (03)(12), crossing (02)(13)}` and let
`ρ_{v,r}(C)` replace the two window lines by local pairing `r`, keeping
`τ`, `q` and every other line.

**Lemma 1 (well-defined re-pairing).** If `W_v` is eligible, the two local
lines hold consecutive first-opening ranks `R, R+1` in every local pairing.
Keeping the `q` tuple fixed therefore assigns `q_R` to the first-opened and
`q_{R+1}` to the second-opened local line. `n`, `τ`, eligibility of `W_v`,
and every line outside `W_v` are invariant under `ρ_{v,r}`, and
`ρ_{v,r'} ∘ ρ_{v,r} = ρ_{v,r'}`.

*Proof.* The window vertices are consecutive positions and all four
belong to local lines. So no other line opens between positions `v` and
`v+3`. Lines opening before `v` have rank `< R` and lines opening after
`v+3` have rank `> R+1`. The other statements follow because `ρ` only
rewrites edges inside `W_v`. ∎

The Jacobian of `ρ_{v,r}` is 1: it is a bijection of discrete labels.
Continuous times and momenta are untouched.

**Lemma 2 (window factorization).** For eligible `W_v`,

```
D(ρ_{v,r} C) = tr[ R_v · W_v^{(r)} · L_v ] · Φ_v,
```

where
- `L_v` is the ordered product up to the segment entering `τ_v`;
- `R_v` is the ordered product from the segment leaving `τ_{v+3}`;
- `Φ_v` holds the phonon factors of the non-local lines;
- `W_v^{(r)} = M_{v+3} U M_{v+2} U M_{v+1} U M_v × (local phonon factors)`,
  evaluated at `k_local = k_v`.

*Proof.* Both local lines close inside the window, so the momentum
entering (`k_v`) equals the momentum leaving it. Spectator lines spanning
the window change only `k_local` and `Φ_v`, and neither depends on `r`.
Matrix order is preserved. ∎

Hence `Σ_r D(ρ_{v,r} C) = tr[R_v Op_v L_v] Φ_v` with `Op_v = Σ_r W_v^{(r)}`.
For `η = 0, ζ = 1`, `Op_v` is **exactly** the compiled Diagram Compiler
two-band `n = 2` operator at `(k_local, q_R, q_{R+1}, τ_{v..v+3})`. This is
where R1 enters the consumer.

## 3. Tiled groups (the partition)

Fix the tiling `v ∈ {0, 4, 8, …}`. Let `E(C)` be the eligible tiled
windows, `m(C) = |E(C)|`, and

```
G(C) = { ρ_{v_1,r_1} ⋯ ρ_{v_m,r_m} C : r ∈ {3 pairings}^m },   |G(C)| = 3^{m(C)}.
```

**Lemma 3 (partition).** `C' ∈ G(C) ⇒ G(C') = G(C)`, so the groups
partition `Ω`.

*Proof.* Tiled windows are disjoint. Eligibility of `W_w` depends only on
the lines touching `W_w`. Re-pairing `W_{w'}` (`w' ≠ w`) rewrites only
lines inside `W_{w'}`, and re-pairing `W_w` keeps it eligible (Lemma 1).
So `E(C') = E(C)`, and the re-pairings of different windows commute. ∎

*Counterexample for a non-tiled rule (negative control NC3).* Take
"the lowest eligible window at any position". Re-pairing `W_v` from nested
to sequential can make `W_{v−2}` eligible, for example when
`(v−2, v−1)` is a line. The "group" then changes under its own move. The
code finds 512 such states at `G = 8, N = 4`.

**Lemma 4 (group factorization).** With `Op_w = Σ_r W_w^{(r)}`,

```
F(G) := Σ_{C∈G} D(C) = tr[ ⋯ Op_{w_2} ⋯ Op_{w_1} ⋯ ] · Φ_G,
```

that is, the chain with every eligible tiled window replaced by its summed
operator. This costs `m` window sums instead of `3^m` diagrams.

*Proof.* Apply Lemma 2 to each tiled window in turn. A window's
`k_local` and the chain between windows do not depend on other windows'
pairings, because each window's local lines close inside it. Linearity
of the trace then gives the product of sums. ∎

## 4. The two measures and their estimators

**A (BOUND_C, native).** `w_A(C) = |Re D(C)|`, `f_O(C) = sgn(Re D(C)) · O(C)`.

**B (SHARED_X_GROUP on tiles).**

```
w_B(C) = |Re F(G(C))| / |G(C)|,
f_O(C) = Σ_{C'∈G(C)} Re D(C') O(C') / |Re F(G(C))|   (= s_G · Ō(G)).
```

**Identity (unbiasedness).** For both measures, `Σ_C w(C) f_O(C) = Σ_C Re D(C) O(C)`.

Hence `⟨O⟩ = E_w[f_O] / E_w[f_1]`.

*Proof for B.* The sum over `C` becomes a sum over groups `G`. Each
group has `|G|` identical terms `|Re F(G)|/|G| · Σ_{C'∈G} Re D(C')O(C')/|Re F(G)|`,
which add up to `Σ_{C∈G} Re D(C) O(C)`. ∎

**Support condition (S).** B needs `Σ_{C∈G} Re D(C) O(C) = 0` whenever
`Re F(G) = 0`. Group-invariant `O` (for example the order `n`) satisfy (S)
automatically. For other `O`, it is a measure-zero condition in
continuous variables, but it must be checked on a finite grid. It holds
in every regime run here: the unbiasedness checks pass to 1e-12.

**Sign theorem.** Define `⟨s⟩_A = Z/Z_A` and `⟨s⟩_B = Z/Z_B`, with
`Z = Σ Re D`, `Z_A = Σ_C |Re D(C)|` and `Z_B = Σ_G |Re F(G)|`.
Since `Re F(G) = Σ_{C∈G} Re D(C)`,

```
Z_B = Σ_G |Σ_{C∈G} Re D(C)| ≤ Σ_G Σ_{C∈G} |Re D(C)| = Z_A,
```

with equality iff every group is sign-coherent. So B **never has a worse
average sign**. It gains exactly where members of a group cancel. For a
positive model (M5's scalar Holstein, or this family at `η = 0`), B has
**no** sign advantage.

Note the doc's option list (§8.2): B samples `|Re F(G)|`, not `|F(G)|`.
Sampling `|F|` would also be legal with the reweighting factor
`Re F/|F|`. Sampling `|D|` while reweighting by `sgn Re D` is **not**
legal (NC2).

## 5. Kernels

Native moves are chosen with state-independent probabilities
`(add/remove 0.30, swap 0.30, time shift 0.25, q change 0.15)`. Each is a
Metropolis–Hastings step on the measure's weight `w`:

| move | forward proposal | reverse |
|---|---|---|
| add | ½ · 1/C(F,2) · 1/L (F free grid points) | remove: ½ · 1/(n+1) |
| remove | ½ · 1/n | add from C' |
| swap (Luo `update_swap`) | 1/(2n−1), adjacent positions `(i, i+1)` exchange line endpoints | same `i`, an involution |
| shift | 1/(2n) · 1/|free points between neighbours| | same set |
| q change | 1/n · 1/(L−1) | symmetric |

Moves that are impossible leave the state unchanged. Group moves use
probability `p_g = 0.3`, replacing native moves at that rate:

* **A heat-bath.** Choose `v` uniformly on `0..2n−4`. This is legal
  because `n` is invariant (Lemma 1). If `W_v` is eligible, draw
  `ρ_{v,r}C` with probability `∝ |Re D|`. This is a Gibbs step on the
  3-element fibre. It needs the three individual `|Re D_r|` values.
  Because the absolute value is not linear, **R1's summed `Op` cannot be
  used here**. Only shared environments `L_v, R_v` (CSE-type reuse) help.
* **B refresh.** Draw `C'` uniformly from `G(C)`. This is Gibbs, because
  `w_B` is constant on groups. It is free: no weight evaluation.

Variable order (§8.2 item 15) needs no group birth/death kernel. The
chain lives on configurations, and add/remove simply recompute `G(C')`
(the tiling is re-derived from the new positions). The multiplicity
`1/|G|` in `w_B` absorbs the change of `m`.

**Negative controls** (each must fail):

| id | construction | expected failure |
|---|---|---|
| NC1 | `w = |Re F(G)|` without `1/|G|` | estimator biased (groups over-weighted by `3^m`) |
| NC2 | sample `|D|`, reweight `sgn Re D` | biased when `D` is complex |
| NC3 | group = first eligible window at any position | not a partition; refresh breaks `πP = π`; estimator biased |
| NC4 | A heat-bath with `v` uniform on **eligible** windows only | the selection probability changes under re-pairing, which breaks `πP = π` |

## 6. What a grouped sampler must define, and where it is defined here

| item | here |
|---|---|
| 1 group definition | §3: tiled eligible windows, all local re-pairings |
| 2 partition vs cover | partition (Lemma 3); cover counterexample NC3 |
| 3 anchor / selection | none needed: B is MH on configurations with a group-constant weight |
| 4 multiplicity | `1/|G| = 3^{−m}` in `w_B`; omitting it is NC1 |
| 5–6 proposal and reverse | §5 table; checked by exact detailed balance |
| 7 target | `|Re F(G)|/|G|`, compared against native `|Re D|` |
| 8 normalization | finite `Ω`, exact enumeration |
| 9 Jacobians | `ρ` is a label bijection (Jacobian 1); grid counting measure |
| 10 sign / phase | complex, k-dependent vertices; sign theorem |
| 11 open legs | spectators spanning a window are carried in `k_local` and `Φ`. External phonon legs are **not modelled** |
| 12 band ordering | all products ordered; non-commuting `U`, `M` |
| 13–14 numerator / denominator | `E_w[f_O]`, `E_w[f_1]` from the same chain |
| 15 variable order | §5 |
| 16 double counting | Lemma 3; NC3 |

## 7. What this does not establish

This work does not establish:
- a LiF calculation;
- native Fortran;
- external phonon legs;
- windows larger than four vertices;
- the continuous-time version, although the proofs carry over with
  Jacobian 1 for `ρ`;
- an evaluator speed comparison;
- equal-precision efficiency.

## 8. Predeclared acceptance (fixed in `run_checks.py` before any run)

* Exact identities and oracle agreement: relative error ≤ `1e-12`.
* Stationarity and detailed balance: `max|πP − π| / max π ≤ 1e-11`;
  support irreducible.
* Every negative control deviates by ≥ `1e-6`, in bias or in stationarity.
* Monte Carlo: 400k steps, 10% burn-in, 50-batch jackknife; each estimate
  within `4σ` of the exact value.
* **Reading of the result.** The construction is *legal* if all of the
  above pass. Whether a material version is *worthwhile* is judged from the
  exact sign gain `Z_A/Z_B` and the exact per-step asymptotic variances.
  These are reported, not tuned.

## 9. Result

`python -m fepdmc_toolkit.r1_grouped.run_checks --out DIR` runs all checks
(117 of 117 pass; needs scipy, which the `dev` and `analysis` extras install). The regression
tests are in `tests/test_r1_grouped.py`.

- The grouped measure is legal: per-diagram weights agree with an
  independent oracle to ≤ 5e-15; the partition (Lemma 3) and factorization
  (Lemma 4) hold, the latter to ≤ 3e-16, and to 1.9e-16 with the compiled
  two-band n = 2 `Op` of the diagram compiler as the window operator; every
  legal kernel is stationary (`max|πP − π| / max π ≤ 2e-15`) and satisfies
  detailed balance to ≤ 2e-17; all negative controls fail as required.
- The sign gain is small: `Z_A/Z_B` is 1.000, 1.073 and 1.068 in the three
  regimes (positive compiler family, signed interband, signed mixed).
- Not carried to a real material in this form; see [roadmap.md](roadmap.md).
