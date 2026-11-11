# Native FEP-DMC run records

JSON summaries of the native runs behind [docs/fepdmc.md](../../docs/fepdmc.md),
written by the analysis modules in `fepdmc_toolkit.fepdmc.analysis` and by
`fepdmc-toolkit compare`. Materials: `lif_elec` (LiF electron), `lif_hole`
(LiF hole), `sto` (SrTiO₃), `anatase`.

| files | contents |
|---|---|
| `external_pair_validation.json` | LiF-hole T = 500 K: maxOrder-7 energies before and after the removal fixes, the published β = 232 cases, full-order sign and pair counts |
| `lif_hole500_allfix_vs_wqfix.json` | full order, all three fixes against the frequency fix only |
| `lif_hole500_ext_pscan.json` | external-move probability scan with all fixes |
| `*_pooled.json`, `*_between.json`, `*_vs_*.json` | pooled comparisons of two run sets (A against B or C): Q difference, variance ratio with 90% interval, CPU or wall ratio, efficiency |
| `medo31_*.json` | maxOrder-31 comparisons against the first 32 native chains. With more batches the net efficiency of general add/remove fell to 1.08 (90% interval 0.67–1.92); a single batch is not enough (see the estimator notes in docs/fepdmc.md) |
| `lif_hole500_long_mix_vs_native.json` | 10⁷-step runs at full order: native against all moves |
| `*_rb*.json`, `*_modes*.json`, `*_mode_rb.json`, `*_svd_basis.json`, `*_qgroup.json`, `*_group_rules_estimate.json`, `*_variance_decomposition.json` | measurement-side diagnostics: window Rao–Blackwellization, phonon-mode statistics, variance decomposition of Q |
| `sto_bchain_*.json` | grouped chain B on SrTiO₃ against native |
| `materials_sign.json`, `materials_sign_lifhole_mo990.json` | average sign of the native sampling per material |
| `BINARY_SHA256.txt` | SHA256 of the `perturbo.x` builds used |
| `FEP_DMC_PIN` | upstream commit of those builds |
