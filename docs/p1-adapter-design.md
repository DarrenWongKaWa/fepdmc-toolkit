# P1 adapter: design and contract

Design of the P1 source adapter in `integration/p1/`. Usage and validation
are in [p1.md](p1.md).

Pinned upstream: `https://github.com/yaoluo/FEP-DMC` commit
`05d08449cffdbd0dfbbbf5009add5cc887bc754b`.
Native subroutine of record: `update_swap` in
`perturbo-fep-dmc/pert-src/diagMC_JJ_updates.f90` (the JJ Method0 path).
Other files named `update_swap` (`diagMC_updates.f90`,
`diagMC_DMFT_updates.f90`, `diagMC_selfe_updates.f90`,
`diagMC_JJt_update.f90`) are out of scope.

The adapter must preserve the native target whatever its wall-clock
gain. The measured gain is in [results.md](results.md).

---

## 1. Scope and non-goals

**In scope**

- Delayed-acceptance P1 on the pinned public JJ `update_swap`.
- Preserve the **native** swap transition kernel’s invariant, as
  implemented by that subroutine.
- Fail closed outside a declared operating domain.
- Compose with the public C0 source adapter without requiring it.
- Reuse the C0 identity pattern: pin → dry-run → exact postimage →
  PRISTINE / APPLIED / UNKNOWN.

**Out of scope**

- P1 on non-JJ `update_swap` copies.
- Neural / linear / LR cheap scores (`LINEAR_DA_SCORE=linear|lr`).
- Claiming ≥5% wall gain or equal-precision efficiency.
- Vendoring FEP-DMC; downloading `lif-sp3_epwan.h5`.
- A generic plugin framework.
- `C2_FORCE_A1`, `P1_FIXTURE`, `c2_dump_swap`, `c1_mark_*` in
  production.

---

## 2. Native public `update_swap` anatomy

Source: pinned `diagMC_JJ_updates.f90`, `subroutine update_swap`
(approximately lines 1015–1163).

Control flow (vanilla):

1. If `diagram%order <= 3`, return.
2. Draw `iv1` uniform in `{2,…,order}`; `iv2 = vL%link(3)` (next
   electron-line vertex).
3. Reject the proposal (return, no RNG accept draw) if:
   - `vL%link(2) == iv2` (phonon partner is that next vertex);
   - `iv2 == maxN`;
   - the pair is the external head/tail crossing pair;
   - `tauL > tauR` (hard stop).
4. Increment `stat%update(7)`.
5. Build trial vertices `vRn` (time of `vL`, phonon of `vR`) and
   `vLn` (time of `vR`, phonon of `vL`):
   - copy `nu`, `i_q`, `Pq`, `wq` from the phonon donor;
   - `vRn%i_kout = vRn%i_kin + vRn%i_q`;
   - `call cal_ek_int` for the new mid-segment electron state;
   - **`call cal_gkq_vtex_int` and `cal_gkq_full_vtex_int`** for both
     trial vertices (expensive).
6. Propagator ratio (cheap, no `g`):
   - `tau12 = tauR - tauL`;
   - `P_kchange = Gel(minval(vRn%ekout), tau12) / Gel(minval(vL%ekout), tau12)`;
   - `wq1 = vL%wq(vL%nu)`, `wq2 = vR%wq(vR%nu)`;
   - multiply by two native `Dph` ratios on the phonon-partner times
     `link(2)`.
7. Environment contraction (expensive):
   - `left_environment_matrix`, `right_environment_matrix`,
     `propagator` / `matmul` for old and new matrix elements;
   - `mat_old = trace(A)` unless `sample_gt`, then
     `A(ib_sample, ib_sample)`.
8. `factor = real(mat_new) / real(mat_old)`
9. `P_accept = abs(factor) * P_kchange`
10. `call random_number_omp(diagram%seed, ran)`
11. If `ran < P_accept`, commit via `copy_vtex` into `iv1`/`iv2` and
    increment `stat%accept(7)`.

**Stage-1 insertion point (required):** after the time-order check
(step 3) and **before** the first `cal_gkq_vtex_int` (step 5). The
benchmarked build places the hook there. A future source test must
fail if the stage-1 `return` occurs after any `cal_gkq_vtex_int` or
`left_environment_matrix` / `right_environment_matrix` in this
subroutine.

**Data available at that point:** existing vertices `vL`, `vR`,
`vLL`, `vRR`; times `tauL`, `tauR`; phonon fields `wq`, `nu`,
`i_q`, `link(2)` partner times; `vL%i_kin`, `vL%ekout`. Not yet
available: trial `gkq`, environment matrices, `mat_old` / `mat_new`.

Stage 1 may call **`cal_ek_int` once** for the proposed mid-segment
`k = vL%i_kin + vR%i_q` (as `p1_prop_ell` in the benchmarked build does). That is
not `g` and not an environment contraction.

`update_swap` copies `wq` from existing vertices. It does not call
`sample_q_omp_int`. Stale `wq` after `add_external_ph` is a C0
defect on a different update.

---

## 3. Mathematical contract

### 3.1 Native ratio

From the pinned body:

```
factor   = real(mat_new) / real(mat_old)
P_accept = abs(factor) * P_kchange
```

Define, for a proposal `x → y` that reaches the accept draw:

\[
R_{\mathrm{native}}(x,y) := P_{\mathrm{accept}}
  = \frac{\lvert \mathrm{Re}\, M_{\mathrm{new}} \rvert}
         {\lvert \mathrm{Re}\, M_{\mathrm{old}} \rvert}\,
    P_{k\mathrm{change}}(x,y)
\]

when \(\mathrm{Re}\, M_{\mathrm{old}}\) is finite and nonzero.
This is **not** \(\lvert M \rvert_{\mathbb{C}}\).

\(P_{k\mathrm{change}}\) is exactly the product in the pinned source:

\[
\begin{aligned}
P_{k\mathrm{change}}
&= \frac{G_{\mathrm{el}}(\min_b e^{\mathrm{new}}_{k,\mathrm{out}},\tau_{12})}
        {G_{\mathrm{el}}(\min_b e^{\mathrm{old}}_{k,\mathrm{out}},\tau_{12})} \\
&\quad\times
  \frac{D_{\mathrm{ph}}(\omega_1,\lvert \tau_R-\tau_{P1}\rvert)}
       {D_{\mathrm{ph}}(\omega_1,\lvert \tau_L-\tau_{P1}\rvert)}
\times
  \frac{D_{\mathrm{ph}}(\omega_2,\lvert \tau_L-\tau_{P2}\rvert)}
       {D_{\mathrm{ph}}(\omega_2,\lvert \tau_R-\tau_{P2}\rvert)}
\end{aligned}
\]

with \(\tau_{12}=\tau_R-\tau_L\), \(\omega_1=v_L\%wq(v_L\%nu)\),
\(\omega_2=v_R\%wq(v_R\%nu)\), \(\tau_{P1}=\) partner time of \(v_L\),
\(\tau_{P2}=\) partner time of \(v_R\).

Pinned primitives (`diagMC.f90`):

- `Gel(ek, tau) = exp(-ek * tau)` (no temperature branch).
- `Dph(w, tau)`:
  - if `w < phfreq_cutoff * ryd2ev`: `1e-15`;
  - else if `zeroTMC`: `exp(-abs(tau) * w)`;
  - else finite-T image sum `D1 + D2` using `config%tauMax`.

`phfreq_cutoff` is read in meV and stored as Rydberg
(`phfreq_cutoff = phfreq_cutoff / ryd2mev`). Default namelist `1.0`
meV therefore matches a `1e-3` eV comparison in `Dph`. A non-default
namelist cutoff must follow native `Dph`, not a hard-coded `1e-3`.

The native kernel for proposals that reach step 10 is
`ran < P_accept`, i.e. Metropolis \(A(x,y)=\min(1,R_{\mathrm{native}})\),
with `ran` from `diagram%seed`. Public P1 does **not** add a Hastings
`q(x,y)/q(y,x)` factor that vanilla omits.

### 3.2 Delayed acceptance

Let \(r_{\mathrm{hat}}(x,y)>0\) be a cheap positive score with

\[
r_{\mathrm{hat}}(y,x) = 1 / r_{\mathrm{hat}}(x,y)
\quad\Leftrightarrow\quad
\ell_{\mathrm{hat}}(y,x) = -\ell_{\mathrm{hat}}(x,y).
\]

\[
\begin{aligned}
\alpha_1(x,y) &= \min\bigl(1, r_{\mathrm{hat}}(x,y)\bigr),
  &\log\alpha_1 &= \min(0,\ell_{\mathrm{hat}}),\\
\alpha_2(x,y) &= \min\bigl(1, R_{\mathrm{native}}(x,y)/r_{\mathrm{hat}}(x,y)\bigr),
  &\log\alpha_2 &= \min(0,\ell_R-\ell_{\mathrm{hat}}).
\end{aligned}
\]

Then \(\alpha_1\alpha_2 / (\alpha_1^{\mathrm{rev}}\alpha_2^{\mathrm{rev}})
= R_{\mathrm{native}}\) whenever the reciprocal-score identity holds
(Christen–Fox). Clipping \(\ell_{\mathrm{hat}}\) to
\([\log(1/10),\log 10]\) preserves reciprocity if the unclipped score
is already antisymmetric, because \(\mathrm{clip}(-\ell)=-\mathrm{clip}(\ell)\).

**Authoritative logs**

- \(\ell_R := \log R_{\mathrm{native}} = \log P_{\mathrm{accept}}\)
  when \(P_{\mathrm{accept}}>0\). Computed from the **native**
  `P_accept` after the vanilla matrix contraction. Do not re-sum
  `Gel`/`Dph` in software for \(\ell_R\).
- \(\ell_{\mathrm{hat}}^{\mathrm{raw}} := \log P_{k\mathrm{change}}\)
  using **native** `Gel` and `Dph` with the same arguments as the
  vanilla `P_kchange` block, except `minval(vRn%ekout)` may be
  obtained from a single `cal_ek_int` without building trial `gkq`.
- \(\ell_{\mathrm{hat}} := \mathrm{clip}(\ell_{\mathrm{hat}}^{\mathrm{raw}},
  \pm\log 10)\).
- \(r_{\mathrm{hat}} = \exp(\ell_{\mathrm{hat}})\).

Stage 1 rejects if
\(\log\max(u_1,10^{-300}) \ge \min(0,\ell_{\mathrm{hat}})\),
with \(u_1\) from the auxiliary stream, and **returns before**
`cal_gkq`.

Stage 2, after native `P_accept` and native `ran`:
accept iff
\(\log\max(\mathrm{ran},10^{-300}) < \min(0,\ell_R-\ell_{\mathrm{hat}})\).

### 3.3 Reverse cheap score

Occupancy reverse at fixed time order: exchange the two phonons’
\((\omega,\tau_{\mathrm{partner}})\); keep \(\tau_L,\tau_R\); exchange
the two mid-segment \(\min e_{k,\mathrm{out}}\) values. Do not flip
both \(\Delta\tau\) and \(\Delta E\).

Under native `Gel` and native `Dph` (including finite-T `D1+D2`,
which depend on `abs(tau)`), that reverse is antisymmetric in
\(\ell_{\mathrm{hat}}^{\mathrm{raw}}\). Clipping preserves it.

---

## 4. Operating domain

Checked once in `linear_da_ensure` after `dmc_band` is known.
`LINEAR_DA=on` with a refused setting is `error stop`, not a silent
native fallback.

| Setting | First public adapter |
|---|---|
| Driver / file | JJ `update_swap` in `diagMC_JJ_updates.f90` only — **SUPPORTED** |
| Other `update_swap` files | **NOT YET IMPLEMENTED** (unpatched; do not claim P1) |
| `DMC_Method` | `0` **SUPPORTED**; any other value **REFUSED** |
| `dmc_band` | `1` **SUPPORTED**; any other value **REFUSED** |
| `sample_gt` | `.false.` **SUPPORTED**; `.true.` **REFUSED** |
| `zeroTMC` | `.true.` **SUPPORTED**; `.false.` **REFUSED** at runtime |
| `LINEAR_DA_SCORE` | `prop` **SUPPORTED**; any other value **REFUSED** |
| `LINEAR_DA` | `off` (default) or `on`; `shadow` **REFUSED** |
| C0 state | independent — see §8 |
| `Holstein` / `LatticeFrohlich` | allowed only if the domain row above still holds after native init |

**Why Method0 and one DMC band:** that is the validated Luo EZ / LiF
setting. DA algebra does not forbid other values if \(\ell_R\) is
native `P_accept`. First release still **refuses** them.

**Why refuse `zeroTMC=.false.` while using native `Dph`:** finite-T
`Dph` in \(\ell_R\) is exact (Architecture C). No public finite-T
target-equivalence fixture exists. Runtime refuse until that fixture
is added; the code path must still call native `Dph` (no hard-coded
`exp(-|τ|ω)` in \(\ell_R\)).

**Why refuse `sample_gt`:** vanilla then replaces `trace(A)` by a
single matrix element. Unvalidated for P1.

---

## 5. Chosen architecture

### Architecture A — transplant of the benchmarked module

Reimplement log `Gel`/`Dph` inside `linear_da_mod` (`linear_ld` with
hard-coded `1e-3` eV and zero-T `−ω|Δτ|`). \(\ell_R\) rebuilt from
those logs plus \(\log\lvert\mathrm{Re}\,M\rvert\).

- Exactness risk: **high** at finite T and non-default
  `phfreq_cutoff`.
- Duplicates the native target definition.
- Carries private symbols.

### Architecture B — independent cheap score, native logs for \(\ell_R\)

Cheap \(\ell_{\mathrm{hat}}\) as in A or as a propagator log;
\(\ell_R\) reconstructed by calling native `Gel`/`Dph` plus
\(\log\lvert\mathrm{Re}\,M\rvert\).

- Exactness: better than A; still duplicates `P_kchange` arithmetic.
- Drift risk if vanilla later changes `P_accept` assembly.

### Architecture C — cheap propagator score; \(\ell_R=\log P_{\mathrm{accept}}\)

Stage 1: \(\ell_{\mathrm{hat}}\) from native `Gel`/`Dph` (log
`P_kchange`) plus one `cal_ek_int`, clipped.
Stage 2: vanilla computes `P_accept` unchanged; \(\ell_R=\log P_{\mathrm{accept}}\).

- Exactness risk: **lowest**. The native target is the native
  variable.
- Duplication: only the cheap score, which is allowed to be
  approximate (clipped).
- Finite T: \(\ell_R\) automatically follows native `Dph`.
- Maintainability: vanilla `factor` / `P_kchange` lines stay the
  single source of truth.
- Performance: same as A/B for the expensive path (stage 1 still
  returns before `cal_gkq`).
- Tests: can assert patched `P_accept = abs(factor)*P_kchange` still
  present, and stage-2 uses `P_accept` not a third formula.

**Chosen: Architecture C.**

Rejected: A (transplant of the benchmarked module). B as
unnecessary duplication of `P_kchange`.

---

## 6. Rejected alternatives (summary)

| Alternative | Why rejected |
|---|---|
| Copy private tree (`c1_mark_*`, `c0_copy_wq`, fixtures) | Vanilla has none of those symbols |
| \(\lvert M\rvert_{\mathbb{C}}\) target | Not the pinned `abs(real(mat))` measure |
| Unclipped cheap score as the accept probability | Changes the target |
| Silent native fallback when domain fails | Hidden retarget / false “P1 on” |
| `C2_FORCE_A1` in production | Breaks \(\ell_{\mathrm{hat}}(y,x)=-\ell_{\mathrm{hat}}(x,y)\) |
| Stage-1 after `cal_gkq` | No skip of expensive work |
| Hard-coded `1e-3` eV in \(\ell_R\) | Conflicts with namelist `phfreq_cutoff` |
| Zero-T `linear_ld` in ell_R | Conflicts with native finite-T Dph |
| Require C0 before P1 | `update_swap` does not sample `q` |
| Generic plugin engine | Out of scope |

---

## 7. Symbols of the benchmarked build

| Symbol in the benchmarked build | Disposition | In the adapter |
|---|---|---|
| `linear_da_is_on` / `linear_da_ensure` | REQUIRED PUBLIC P1 CORE | Keep, domain-checked |
| `linear_da_eval` + `p1_prop_ell` (prop score) | REQUIRED PUBLIC P1 CORE | Rewrite: native `Gel`/`Dph` logs, no `1e-3` hardcode in \(\ell_R\) |
| `linear_da_aux_uniform` | REQUIRED PUBLIC P1 CORE | Keep independent xorshift; default seed `9142871` |
| clip \(\pm\log 10\) | REQUIRED PUBLIC P1 CORE | Keep; not a namelist |
| `n_invoked`, `n_eligible`, `n_stage1_rej`, `n_stage2`, `n_accept` | OPTIONAL OBSERVABILITY | Keep a reduced set (§14); print `LINEAR_DA_COUNTS` |
| `linear_da_report` | OPTIONAL OBSERVABILITY | Keep |
| `LINEAR_DA=shadow` | REMOVE | Refuse |
| `LINEAR_DA_SCORE=linear\|lr`, `prefix_logp`, `ids_closed`, `ngrid_da` | REMOVE | Refuse non-`prop` |
| `p1_fixture_on` | TEST ONLY / REMOVE from production | Python tests; no Fortran fixture gate |
| `C2_FORCE_A1` | REMOVE | Forbidden in production source |
| `c2_dump_swap`, `c2_live_fp`, `C2_DUMP` | BENCHMARK BUILD ONLY | Remove |
| `c1_mark_*`, `c1_ncall`, `closure_observer` | BENCHMARK BUILD ONLY | Remove; counters in `linear_da_mod` |
| `c0_copy_wq` | BENCHMARK BUILD ONLY | Native `vRn%wq = vR%wq` assignment |
| `src/fepdmc_toolkit/fortran/linear_da_mod.f90` | reference copy of the benchmarked module | Do not copy into vanilla; the adapter's module is `integration/p1/linear_da_mod.f90` |
| `src/fepdmc_toolkit/fortran/update_swap_p1_excerpt.f90` | reference copy of the benchmarked excerpt | Reference only |
| `fepdmc_toolkit.p1_da` | Python algebra used by the tests | Not linked into `perturbo.x` |

---

## 8. C0 composition

**Policy: independent, composable, either order.**

`update_swap` never calls `sample_q_omp_int`. C0 repairs
`add_external_ph` in the **same file**, different subroutine.
P1 matching native swap does not require C0. A chain that used
`add_external_ph` without C0 already has native stale-`wq` error;
P1 does not add a second error source.

C0 and P1 both edit
`perturbo-fep-dmc/pert-src/diagMC_JJ_updates.f90`, so whole-file
hashes are not independent. Record **four** full-file SHA256 values
for the pinned commit (recorded in `provenance/P1_PUBLIC_PATCH.json`):

| State name | Meaning |
|---|---|
| `PUBLIC_PRISTINE` | vanilla pin preimage |
| `PUBLIC_C0_APPLIED` | C0 postimage |
| `PUBLIC_P1_APPLIED` | P1 applied to vanilla |
| `PUBLIC_C0_P1_APPLIED` | both, either order |
| `UNKNOWN` | anything else |

P1 apply:

- `PUBLIC_PRISTINE` → transform P1 region → `PUBLIC_P1_APPLIED`
- `PUBLIC_C0_APPLIED` → same P1 region transform → `PUBLIC_C0_P1_APPLIED`
- `PUBLIC_P1_APPLIED` or `PUBLIC_C0_P1_APPLIED` → `ALREADY_APPLIED`, exit 2
- else refuse

C0 apply:

- `PUBLIC_PRISTINE` → C0 (already specified)
- `PUBLIC_P1_APPLIED` → C0 region transform → `PUBLIC_C0_P1_APPLIED`
- `PUBLIC_C0_APPLIED` or `PUBLIC_C0_P1_APPLIED` → already applied
- else refuse

P1 must not rewrite the C0 insertion. C0 must not rewrite the P1
insertion. Anchors live in different subroutines; tests must fail if
either transform touches the other subroutine’s body.

---

## 9. Source transformation design

Files:

```
integration/p1/apply_p1.py
integration/p1/verify_p1.py
integration/p1/p1_lib.py
integration/p1/linear_da_mod.f90
integration/p1/README.md
provenance/P1_PUBLIC_PATCH.json
```

**Preimage states accepted by apply:** `PUBLIC_PRISTINE` and
`PUBLIC_C0_APPLIED` as defined by
`provenance/UPSTREAM_FEP_DMC.json` and
`provenance/C0_PUBLIC_PATCH.json`.

**Transforms (mechanical)**

1. Copy `integration/p1/linear_da_mod.f90` to
   `perturbo-fep-dmc/pert-src/linear_da_mod.f90` if absent or if the
   bytes equal the canonical module; refuse if a different file
   occupies that path.
2. Makefile: in `perturbo-fep-dmc/pert-src/makefile`, unique
   `PERTMOD` fragment

   ```
   diagMC.f90 \
   diagMC_debug.f90 \
   ```

   becomes

   ```
   diagMC.f90 \
   linear_da_mod.f90 \
   diagMC_debug.f90 \
   ```

   Refuse if the fragment is missing or duplicated. `OBJ` is derived
   from `PERTMOD`; no other makefile line changes.
3. `update_swap` in `diagMC_JJ_updates.f90` only:
   - add `use linear_da_mod, only:` the production symbols in §10;
   - after the `tauL>tauR` stop, insert Stage 1 (eval, aux uniform,
     possible `return` before any `cal_gkq`);
   - keep vanilla trial construction and `P_accept` assembly;
   - after `random_number_omp`, if P1 on: Stage 2 from `P_accept`;
     if P1 off: `da_acc = (ran < P_accept)`;
   - commit iff `da_acc`.
4. `diagMC_JJ.f90` as specified in section 14: `use linear_da_mod`,
   `call linear_da_ensure()` after `setup_dqmc()`, 
   `call linear_da_report()` after the acceptance write.

Dry-run: same validation, print unified diffs of makefile +
`diagMC_JJ_updates.f90` + `diagMC_JJ.f90` + “would add
`linear_da_mod.f90`”, zero writes.

Atomic writes: temp file in the same directory, `os.replace`,
preserve mode (C0 lesson: public Fortran is `100755`).

Rollback: `git restore` on the three modified files and
`rm perturbo-fep-dmc/pert-src/linear_da_mod.f90` if it was added.
Do not run git restore automatically.

Idempotence: second apply on an exact P1 postimage of either C0
combo exits 2 without writes.

Upstream verifier (`tools/verify_fep_dmc_upstream.py`) remains
pristine-only and **must fail** on any P1-applied tree.

---

## 10. Fortran module API design

Production module name: `linear_da_mod`
(file `linear_da_mod.f90` in `pert-src` after apply).

`use DiagMC` (types, `Gel`, `Dph`, `cal_ek_int`, `dmc_band`).
`use pert_param` (`DMC_Method`, `zeroTMC`, `phfreq_cutoff` only as
consumed inside native `Dph`, not reimplemented).
`use pert_const` (`dp`).

Public procedures:

```
subroutine linear_da_ensure()
  ! Read LINEAR_DA, LINEAR_DA_AUX_SEED.
  ! error stop if LINEAR_DA is not off/on.
  ! error stop if LINEAR_DA=on and domain row in §4 fails.
  ! Initialize aux RNG. Idempotent.

logical function linear_da_is_on()
  ! .true. iff LINEAR_DA=on after ensure.

subroutine linear_da_eval(diagram, iv1, iv2, eligible, ell_hat)
  ! prop score only. eligible is .true. whenever domain passed
  ! and order>3 (the caller already enforced order and topology).
  ! ell_hat = clip(log P_kchange) via native Gel/Dph + one cal_ek_int.
  ! Increments n_seen and n_eligible.

subroutine linear_da_aux_uniform(u)
  ! xorshift64; does not touch diagram%seed.

subroutine linear_da_stage2(ran, P_accept, ell_hat, da_acc)
  ! ell_R = log(P_accept) if P_accept > 0.
  ! Policy §13 for zeros / NaN / Inf.
  ! da_acc = log(max(ran,1e-300)) < min(0, ell_R - ell_hat)

subroutine linear_da_note_expensive()
  ! n_expensive = n_expensive + 1; called immediately before
  ! the first cal_gkq_vtex_int in the patched update_swap.

subroutine linear_da_note_accept()
subroutine linear_da_note_stage2_reject()
subroutine linear_da_report()
  ! prints LINEAR_DA_COUNTS ...
```

No `p1_fixture_on`, no `C2_FORCE_A1`, no dump TSV, no
`closure_observer`.

`src/fepdmc_toolkit/fortran/linear_da_mod.f90` is the reference copy of
the benchmarked module and is **not** the file copied into vanilla.

---

## 11. Configuration

| Name | Kind | Default | Legal values |
|---|---|---|---|
| `LINEAR_DA` | environment | `off` | `off`, `on` |
| `LINEAR_DA_AUX_SEED` | environment | `9142871` | nonzero 64-bit integer; `0` becomes `1` |
| `LINEAR_DA_SCORE` | environment | `prop` | `prop` only |

No new `perturbo` namelist keys in the first public adapter (avoids
editing `pert_param.f90`). Domain quantities `DMC_Method`,
`dmc_band`, `zeroTMC`, `sample_gt` are native namelist/init fields
already.

Missing `LINEAR_DA` → `off`.
Unsupported value → `error stop` at `linear_da_ensure`.
`LINEAR_DA=on` is validated in `linear_da_ensure` called from
`diagmc_JJ` after `setup_dqmc` (section 14). Do not initialize from
a static `block data` unit.

---

## 12. RNG policy

| Draw | Generator | When |
|---|---|---|
| Stage 1 \(u_1\) | auxiliary xorshift64 (`linear_da_aux_uniform`) | P1 on and eligible, before `cal_gkq` |
| Stage 2 `ran` | native `random_number_omp(diagram%seed, ran)` | every swap that reaches the vanilla accept line |

**P1 OFF guarantee (precise):** for every invocation of JJ
`update_swap`, `linear_da_is_on()` is false. The `if (linear_da_is_on())`
guards around Stage 1 and Stage 2 do not run their bodies. No
auxiliary draw. No extra `cal_ek_int` from P1. Native
`diagram%seed` consumption on that subroutine equals vanilla. The
accept test is `ran < P_accept`. This is **same transition kernel
and same native RNG sequence** for that update, not a claim of
whole-binary bit-identity across compilers.

**P1 ON:** stage-1 rejects do **not** consume `diagram%seed`. Vanilla
would have consumed one `ran` after expensive work. Different native
stream; **same target** via DA, not the same sample path.

Aux seed default `9142871` matches the benchmarked module so
Python algebra tests that hard-code that seed remain meaningful.
Changing the aux seed changes P1-on samples, not P1-off samples.

---

## 13. Numerical / error policy

| Event | Policy |
|---|---|
| Unsupported domain + `LINEAR_DA=on` | `error stop` in `ensure` |
| `Re(mat_old)` zero or non-finite | `error stop` with a fixed message; do not accept `Inf` as native `ran < Inf` would | 
| `Re(mat_new) = 0` with finite nonzero old | `da_acc = .false.` (native `P_accept=0`) |
| `P_accept <= 0` | stage-2 reject; do not take `log` of a non-positive number |
| `P_accept` non-finite | `error stop` |
| `u_1` or `ran` non-finite or `<= 0` | floor to `1e-300` for the log test only (representation safeguard, same as `max(ran,1e-300)` in the benchmarked build); this is not a change of \(R_{\mathrm{native}}\) |
| `Dph` returns `1e-15` (native cutoff) | cheap log uses `log(1e-15)`; exact path uses native `P_kchange` |
| Underflow in `Gel`/`Dph` products | native `P_accept` may underflow to 0 → stage-2 reject |
| Overflow in `P_accept` | `error stop` if non-finite; if finite `>1`, native and DA both accept with probability 1 at stage 2 when \(\ell_R-\ell_{\mathrm{hat}}\ge 0\) |
| \(\ell_{\mathrm{hat}}\) outside \(\pm\log 10\) before clip | clip for stage 1 only; never clip \(\ell_R\) |
| NaN in cheap `cal_ek_int` energies | `error stop` |

The `Re(mat_old)=0` stop **differs** from vanilla (`factor` would be
`Inf` and `ran < Inf` accepts). That set is treated as a numerical
singularity, not as a physical accept. Document it in the
implementation README. Do not insert an arbitrary positive floor
into \(R_{\mathrm{native}}\).

---

## 14. Observability

Integer counters in `linear_da_mod`, incremented even when P1 is on
and a path returns early. They do not enter `P_accept`.

| Counter | Increments when |
|---|---|
| `n_seen` | `linear_da_eval` called |
| `n_eligible` | domain-ok swap, order>3, topology already passed |
| `n_s1_reject` | stage-1 return |
| `n_expensive` | about to execute the first `cal_gkq_vtex_int` |
| `n_s2_reject` | stage-2 reject after native `P_accept` |
| `n_accept` | commit |

Invariant to test when P1 is on and the process exits cleanly:

`n_s1_reject + n_expensive = n_eligible`

and `n_expensive = n_accept + n_s2_reject`.

That is the skip proof: stage-1 rejects never call `cal_gkq`.

`diagMC_JJ.f90` (pinned `subroutine diagmc_JJ`) is the report and
init site. Two insertions, no other driver file:

1. After `call setup_dqmc()` (pinned line 26), before the OpenMP
   parallel region: `call linear_da_ensure()`. `dmc_band` is already
   set inside `setup_dqmc`. This avoids first-touch races on
   `da_cfg_init`.
2. After
   `write(stdout,'(A30,7E15.5)')'acceptance = ',...`
   (pinned line 98) and before `call stop_clock('diagmc_JJ')`
   (pinned line 99): `call linear_da_report()`. Same stdout pattern
   as the existing acceptance line (no extra `ionode` guard).

Add `use linear_da_mod, only: linear_da_ensure, linear_da_report`
to `diagmc_JJ`.

OpenMP: JJ `update_swap` runs inside `!$omp parallel`. Module
counters use `!$omp atomic update` on every increment. `ensure` is
not called from the parallel region.

No TSV, no `c1_mark_*`. Profiling clocks are out of scope for the
first public adapter.

---

## 15. Build integration

Pinned `perturbo-fep-dmc/pert-src/makefile`:

- `PERTMOD` compiles as modules in listed order.
- `linear_da_mod.f90` `use`s `DiagMC` (`diagMC.f90`) and
  `pert_param` / `pert_const` (already earlier in `PERTMOD`).
- `diagMC_JJ_updates.f90` will `use linear_da_mod`.
- Insert `linear_da_mod.f90 \` immediately after `diagMC.f90 \`.
- `OBJ = $(PERTMOD:.f90=.o) $(PERTSRC:.f90=.o)` picks up the new
  `.o` without a second edit.
- Link line unchanged.

Compiler: whatever the public `make.sys` already uses. No new
flags. Intel vs gfortran is an environment issue, not an adapter
issue.

---

## 16. Source / provenance states

See §8. Machine record `provenance/P1_PUBLIC_PATCH.json`, fields:

- `schema_version`, `adapter=p1_delayed_acceptance`
- `upstream_commit` (40-char pin)
- `architecture=C`
- `preimage_file` for JJ updates and makefile
- `module_file=perturbo-fep-dmc/pert-src/linear_da_mod.f90`
- `module_sha256` (canonical bytes)
- `states.PUBLIC_PRISTINE`, `PUBLIC_C0_APPLIED`,
  `PUBLIC_P1_APPLIED`, `PUBLIC_C0_P1_APPLIED` each listing SHA256
  for `diagMC_JJ_updates.f90`, `makefile`, and `diagMC_JJ.f90`
- `historical_donor_equivalence=not_established`
- `build_level_established` filled by implementation evidence

Hashes are computed from a live checkout of the pin.

---

## 17. Testing ladder

| Level | Proves | Does not prove | Dependencies |
|---|---|---|---|
| **0** source transform | exact files/hashes; Stage 1 before `cal_gkq`; makefile insert; idempotence; UNKNOWN refuse; C0 compose; no `C2_FORCE_A1` in production module text; `P_accept = abs(factor)*P_kchange` still present | compile, runtime, target | Future B + pin clone |
| **1** compile `linear_da_mod.o` + patched `diagMC_JJ_updates.o` | modules resolve | full `perturbo.x`, physics | `gfortran`/`ifort`, QE modules as public build already needs |
| **2** `perturbo.x` link | the binary contains the module | correct kernel | full public Perturbo+QE build |
| **3A** tiny native fixture (no LiF HDF5): one swap toy or Holstein if public examples allow **without** private paths | P1 off matches vanilla accept bit on that fixture; P1 on stage-1 skip count; domain refuse | LiF Q, HAC, wall ratio | public `example/` only if it builds; otherwise a synthetic `Dph`/`Gel` unit in Fortran tests under Future B that **does not** claim to be `perturbo.x` |
| **3B** public LiF runtime | a real material run can start | statistical exactness | `lif-sp3_epwan.h5` (not included); **not** required to ship P1 source adapter |
| **4** target-equivalence experiment | two-sample test that P1-on and P1-off estimates of a declared observable agree within a **preregistered** protocol | wall-clock gain; six-chain Q agreement alone is **not** this test | LEVEL 2 binary + declared observable |

**Falsification tests (LEVEL 0 / Python):**

- replace `abs(real(M))` by `abs(M)` → algebra test fails
  (`tests/test_da_balance.py` already encodes this; keep it);
- drop `P_kchange` from \(R_{\mathrm{native}}\) → stage-2 identity
  test fails;
- drop `− ell_hat` in stage 2 → `test_p1_clean_source` equivalent
  on the **patched** public file fails;
- `zeroTMC=.false.` with `LINEAR_DA=on` → `error stop` / adapter
  unit test of `ensure`;
- `DMC_Method/=0` or `dmc_band/=1` → refuse;
- wrong upstream commit → apply refuses;
- makefile fragment missing → apply refuses;
- `C2_FORCE_A1` string in production `integration/p1/linear_da_mod.f90`
  → source gate fails;
- stage-1 `return` after `cal_gkq_vtex_int` → source gate fails;
- C0 body bytes changed by P1 transform → compose test fails.

**DA-off equivalence test:** LEVEL 3A (or LEVEL 2 with a public
example) run of identical `diagMC.in` / seed with `LINEAR_DA=off`
on the P1-patched binary versus vanilla binary: native accept/reject
sequence on `update_swap` matches. If a full binary pair cannot be
built, LEVEL 0 proves source-equivalence of the off path (the
vanilla accept line is unchanged and not preceded by aux RNG).

**DA-on target test:** LEVEL 4 protocol, not six-chain Q. Until
LEVEL 4 runs, the implementation may only claim DA algebra + native
`P_accept` wiring.

---

## 18. Remaining scientific risks

These remain after this design. They are facts about evidence, not
holes in the adapter contract.

1. Vanilla Metropolis omits an explicit Hastings ratio. P1 matches
   vanilla; it does not prove vanilla is the diagram-weight
   MH kernel.
2. `Re(mat_old)=0` policy differs from vanilla `Inf` accept. Expected
   measure is negligible; unproven on LiF traces.
3. Finite-T exactness is architectural (native `Dph` in
   `P_accept`) but **runtime-refused** until a finite-T fixture
   exists.
4. Multiband / `DMC_Method/=0` / `sample_gt` unvalidated.
5. The pinned public tree is not proven identical to the tree of the
   benchmarked build (`historical_donor_equivalence=not_established`),
   so an adapter build is not expected to reproduce the benchmark's
   wall-time ratio.
6. LEVEL 3B/4 need data and a compiler stack this repository does
   not ship.
7. OpenMP atomic counters are required; a non-atomic implementation
   makes the skip invariant untestable under the JJ parallel loop.
