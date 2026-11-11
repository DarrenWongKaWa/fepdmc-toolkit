"""Raw vs window-Rao-Blackwellized Luo EZ energy estimator on native chains.

Reads <runs>/<material>/chain*/rb_trace.dat (one line per measurement: raw and
RB numerator/denominator), plus rb_summary.dat self-check counters and the
native parameter.dat-1 energy. For each chain:

    E = sum(num) / (tau_max * sum(den)),  Q = E - E_bare

with a 50-block jackknife standard error for raw and RB on the same
measurements (paired). Reports per-chain SE ratios and the pooled variance
ratio. Predeclared reading: RB is useful only if the pooled SE^2 ratio is
clearly below 1 and the RB energy agrees with raw within errors.

Traces from the general rb_window.f90 also carry rho_single and rho_tiled.
Under the native measure E[rho_tiled] = Z_B/Z_A, so the grouped measure B has
<s>_B = <s>_A / E[rho_tiled] and its sign-limited variance gain is 1/E[rho]^2.
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np

from fepdmc_toolkit.formation_energy import E_BARE_LIF_METHOD0

from ._common import TAU_MAX_50K, AnalysisError, emit, make_parser

NBLOCK = 50


def jackknife_ratio(num, den, tau_max, nblock=NBLOCK):
    nb = len(num) // nblock
    bn = num[: nb * nblock].reshape(nblock, nb).sum(axis=1)
    bd = den[: nb * nblock].reshape(nblock, nb).sum(axis=1)
    est = bn.sum() / bd.sum() / tau_max
    jk = np.array([(bn.sum() - bn[i]) / (bd.sum() - bd[i]) / tau_max for i in range(nblock)])
    se = math.sqrt((nblock - 1) / nblock * np.sum((jk - jk.mean()) ** 2))
    return float(est), float(se)


def read_native_energy(chain: Path):
    p = chain / "parameter.dat-1"
    if not p.exists():
        return None
    lines = p.read_text().split()
    for i, tok in enumerate(lines):
        if tok == "#Etrue":
            return float(lines[i + 2])
    return None


def analyze(runs: str, material: str, tau_max: float = TAU_MAX_50K, e_bare: float = E_BARE_LIF_METHOD0) -> dict:
    root = Path(runs) / material
    rows = []
    for chain in sorted(root.glob("chain*")):
        tr = chain / "rb_trace.dat"
        if not tr.exists():
            continue
        a = np.loadtxt(tr)
        summ = dict(l.split() for l in (chain / "rb_summary.dat").read_text().splitlines())
        extra = {}
        if a.shape[1] >= 10:
            rs, rt = a[:, 7], a[:, 8]
            extra = dict(
                E_rho_single=float(rs.mean()),
                E_rho_tiled=float(rt.mean()),
                frac_meas_rho_tiled_lt_1=float(np.mean(rt < 1 - 1e-12)),
                mean_m_tiled=float(a[:, 9].mean()),
                frac_capped=float(np.mean(a[:, 9] > a[:, 10])) if a.shape[1] >= 11 else None,
            )
        e_raw, se_raw = jackknife_ratio(a[:, 0], a[:, 1], tau_max)
        e_rb, se_rb = jackknife_ratio(a[:, 2], a[:, 3], tau_max)
        rows.append(
            dict(
                chain=chain.name,
                n_meas=len(a),
                E_raw=e_raw,
                se_raw=se_raw,
                E_rb=e_rb,
                se_rb=se_rb,
                se2_ratio=(se_rb / se_raw) ** 2,
                E_native=read_native_energy(chain),
                mean_sign_raw=float(a[:, 1].mean()),
                mean_sign_rb=float(a[:, 3].mean()),
                eligible_fraction=float(a[:, 4].sum() / max(a[:, 5].sum(), 1)),
                mean_order=float(a[:, 6].mean()),
                **extra,
                selfcheck={
                    k: v
                    for k, v in summ.items()
                    if "mismatch" in k or "rel" in k or k in ("n_zero_fibre", "n_tiled_capped")
                },
            )
        )
    if not rows:
        raise AnalysisError(f"{root}: no chain with rb_trace.dat")
    ratios = np.array([r["se2_ratio"] for r in rows])
    Er = np.array([r["E_raw"] for r in rows])
    Eb = np.array([r["E_rb"] for r in rows])
    n = len(rows)
    out = {
        "material": material,
        "tau_max": tau_max,
        "chains": n,
        "per_chain": rows,
        "pooled": {
            "E_bare": e_bare,
            "Q_raw": float(Er.mean() - e_bare),
            "se_Q_raw_between_chains": float(Er.std(ddof=1) / math.sqrt(n)),
            "Q_rb": float(Eb.mean() - e_bare),
            "se_Q_rb_between_chains": float(Eb.std(ddof=1) / math.sqrt(n)),
            "mean_se2_ratio_rb_over_raw": float(ratios.mean()),
            "se_of_mean_se2_ratio": float(ratios.std(ddof=1) / math.sqrt(n)),
            "between_chain_var_ratio": float(Eb.var(ddof=1) / Er.var(ddof=1)),
            "max_abs_E_rb_minus_raw": float(np.max(np.abs(Eb - Er))),
            "native_matches_raw": bool(
                all(r["E_native"] is None or abs(r["E_native"] - r["E_raw"]) < 1e-6 * abs(r["E_raw"]) for r in rows)
            ),
            "selfchecks_clean": all(
                int(v) == 0 for r in rows for k, v in r["selfcheck"].items() if k.endswith("mismatch")
            ),
        },
    }
    if "E_rho_tiled" in rows[0]:
        sA = np.array([r["mean_sign_raw"] for r in rows])
        rt = np.array([r["E_rho_tiled"] for r in rows])
        rs = np.array([r["E_rho_single"] for r in rows])
        out["pooled"].update(
            {
                "sign_A": float(sA.mean()),
                "se_sign_A": float(sA.std(ddof=1) / math.sqrt(n)),
                "E_rho_tiled": float(rt.mean()),
                "se_E_rho_tiled": float(rt.std(ddof=1) / math.sqrt(n)),
                "sign_B_tiled": float(sA.mean() / rt.mean()),
                "sign_gain_B_over_A": float(1.0 / rt.mean()),
                "variance_gain_from_sign_B": float(1.0 / rt.mean() ** 2),
                "E_rho_single": float(rs.mean()),
                "variance_gain_from_sign_single_window": float(1.0 / rs.mean() ** 2),
                "variance_ceiling_perfect_grouping": float(1.0 / sA.mean() ** 2),
                "fraction_of_sign_deficit_recovered": float((sA.mean() / rt.mean() - sA.mean()) / (1 - sA.mean())),
            }
        )
    return out


def main(argv: list[str] | None = None) -> int:
    ap = make_parser(__name__, __doc__)
    ap.add_argument("runs", help="runs with rb_trace.dat (window Rao-Blackwell diagnostic)")
    ap.add_argument("material")
    ap.add_argument("tau_max", type=float, nargs="?", default=TAU_MAX_50K, help="1/kT in 1/eV (default: 50 K)")
    ap.add_argument(
        "e_bare",
        type=float,
        nargs="?",
        default=E_BARE_LIF_METHOD0,
        help="bare band energy in eV, the native 'Ek-mu' (default: LiF electron)",
    )
    a = ap.parse_args(argv)
    return emit(lambda: analyze(a.runs, a.material, a.tau_max, a.e_bare))


if __name__ == "__main__":
    sys.exit(main())
