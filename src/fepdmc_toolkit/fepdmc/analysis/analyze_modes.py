"""Phonon-mode grouping diagnostics on native chains (mode_trace.dat + rb_trace.dat).

Exact under the native measure (E[rho] = Z_B/Z_A for each nu-independent partition):
  single line (mean over lines), first 2 / 3 / 4 lines jointly.
Independence test: E[rho_joint(K)] vs E[prod of the K single-line rho].
Estimate (validated only as far as the test holds): E[prod over all lines of the
single-line rho] ~ Z_B/Z_A for the group that sums the modes of every line.
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np

from ._common import AnalysisError, emit, make_parser


def mean_se(x):
    x = np.asarray(x, float)
    return float(x.mean()), float(x.std(ddof=1) / math.sqrt(len(x))) if len(x) > 1 else float("nan")


def analyze(runs: str, material: str) -> dict:
    rows = []
    for chain in sorted((Path(runs) / material).glob("chain*")):
        mt, rt = chain / "mode_trace.dat", chain / "rb_trace.dat"
        if not (mt.exists() and rt.exists()):
            continue
        m = np.loadtxt(mt)
        r = np.loadtxt(rt)
        n = min(len(m), len(r))
        m, r = m[:n], r[:n]
        sA = r[:, 1].mean()
        summ = dict(l.split() for l in (chain / "rb_summary.dat").read_text().splitlines())
        row = dict(
            chain=chain.name,
            n_meas=n,
            sign_A=float(sA),
            mean_lines=float(m[:, 0].mean()),
            E_rho_tiled4=float(r[:, 8].mean()),
            E_rho_line=float(m[:, 1].mean()),
            E_rho_K2=float(m[:, 3].mean()),
            E_rho_K3=float(m[:, 4].mean()),
            frac_lines_lt_099=float(m[:, 5].mean()),
            selfcheck={k: v for k, v in summ.items() if "mismatch" in k},
        )
        if m.shape[1] >= 11:
            l1, l2, l3 = m[:, 6], m[:, 7], m[:, 8]
            row.update(
                E_rho_l1=float(l1.mean()),
                E_rho_l2=float(l2.mean()),
                E_rho_l3=float(l3.mean()),
                E_prod_l1l2=float((l1 * l2).mean()),
                E_prod_l1l2l3=float((l1 * l2 * l3).mean()),
                E_rho_K4=float(m[m[:, 9] >= 0, 9].mean()) if (m[:, 9] >= 0).any() else None,
                n_K4=int((m[:, 9] >= 0).sum()),
                E_prod_all=float(m[:, 10].mean()),
            )
        if m.shape[1] >= 17:
            lam_n, lam_p, rel, pi, pj, pij = m[:, 11], m[:, 12], m[:, 13], m[:, 14], m[:, 15], m[:, 16]
            row.update(
                E_prod_laminar=float(lam_p.mean()),
                laminar_fraction_of_lines=float(np.mean(lam_n / np.maximum(m[:, 0], 1))),
            )
            for code, name in ((0, "disjoint"), (1, "nested"), (2, "crossing")):
                sel = rel == code
                if sel.any():
                    row[f"pair_{name}"] = dict(
                        n=int(sel.sum()),
                        E_joint=float(pij[sel].mean()),
                        E_product=float((pi[sel] * pj[sel]).mean()),
                        joint_over_product=float(pij[sel].mean() / (pi[sel] * pj[sel]).mean()),
                        E_deficit_joint=float(1 - pij[sel].mean()),
                        E_deficit_sum=float((1 - pi[sel]).mean() + (1 - pj[sel]).mean()),
                    )
        rows.append(row)
    if not rows:
        raise AnalysisError(f"{Path(runs) / material}: no chain with mode_trace.dat and rb_trace.dat")
    out = {"material": material, "chains": len(rows), "per_chain": rows, "pooled": {}}
    keys = [k for k in rows[0] if (k.startswith("E_") or k in ("sign_A", "mean_lines")) and rows[0][k] is not None]
    for k in keys:
        mu, se = mean_se([r[k] for r in rows])
        out["pooled"][k] = {"mean": mu, "se": se}
    p = out["pooled"]
    sA = p["sign_A"]["mean"]
    derived = {}
    for name in ("E_rho_tiled4", "E_rho_line", "E_rho_K2", "E_rho_K3", "E_rho_K4", "E_prod_all", "E_prod_laminar"):
        if name in p:
            rho = p[name]["mean"]
            derived[name] = {
                "sign_B": sA / rho,
                "variance_gain": 1 / rho**2,
                "fraction_of_sign_deficit_recovered": (sA / rho - sA) / (1 - sA),
            }
    derived["ceiling_variance_gain"] = 1 / sA**2
    if "E_prod_l1l2" in p:
        derived["independence_K2_joint_over_product"] = p["E_rho_K2"]["mean"] / p["E_prod_l1l2"]["mean"]
        derived["independence_K3_joint_over_product"] = p["E_rho_K3"]["mean"] / p["E_prod_l1l2l3"]["mean"]
    out["derived"] = derived
    return out


def main(argv: list[str] | None = None) -> int:
    ap = make_parser(__name__, __doc__)
    ap.add_argument("runs", help="runs made with --modes (mode_trace.dat, rb_trace.dat)")
    ap.add_argument("material")
    a = ap.parse_args(argv)
    return emit(lambda: analyze(a.runs, a.material))


if __name__ == "__main__":
    sys.exit(main())
