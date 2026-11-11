"""Raw native estimator vs measurement-side mode Rao-Blackwell on the same chains.

mode_rb_trace.dat: raw_num raw_den mrb_num mrb_den n_lines (one row per
measurement). Both estimators come from the same native trajectory, so the
comparison is paired: per-chain 50-block jackknife SE^2 ratio, and the
between-chain variance of Q. Cost: CPU user time vs hook-free native chains
with the same seeds (identical trajectories).
"""

import math
import re
import sys
from pathlib import Path

import numpy as np

from ._common import TAU_MAX_50K, AnalysisError, emit, make_parser

NB = 50


def jk(num, den, tau):
    nb = len(num) // NB
    bn = num[: nb * NB].reshape(NB, nb).sum(1)
    bd = den[: nb * NB].reshape(NB, nb).sum(1)
    e = bn.sum() / bd.sum() / tau
    j = np.array([(bn.sum() - bn[i]) / (bd.sum() - bd[i]) / tau for i in range(NB)])
    return float(e), float(math.sqrt((NB - 1) / NB * np.sum((j - j.mean()) ** 2)))


def user(c):
    t = re.search(r"user\s+(\d+)m([\d.]+)s", (c / "stderr.log").read_text())
    return int(t.group(1)) * 60 + float(t.group(2)) if t else float("nan")


def analyze(runs: str, plain: str, mat: str, eb: float, tau: float = TAU_MAX_50K) -> dict:
    from scipy import stats  # optional dependency: pip install 'fepdmc-toolkit[analysis]'

    rows = []
    for c in sorted((Path(runs) / mat).glob("chain*")):
        f = c / "mode_rb_trace.dat"
        if not (c / "parameter.dat-1").exists():
            continue
        a = np.loadtxt(f)
        er, ser = jk(a[:, 0], a[:, 1], tau)
        em, sem = jk(a[:, 2], a[:, 3], tau)
        summ = dict(l.split(None, 1) for l in (c / "rb_summary.dat").read_text().splitlines())
        pc = Path(plain) / mat / c.name
        rows.append(
            dict(
                chain=c.name,
                Q_raw=er - eb,
                se_raw=ser,
                Q_mrb=em - eb,
                se_mrb=sem,
                se2_ratio=(ser / sem) ** 2,
                sign_raw=float(a[:, 1].mean()),
                sign_mrb=float(a[:, 3].mean()),
                user_s=user(c),
                user_plain_s=user(pc) if pc.exists() else float("nan"),
                mode_rb_mismatch=summ.get("n_mode_rb_mismatch", "?").strip(),
            )
        )
    if len(rows) < 2:
        raise AnalysisError(f"{Path(runs) / mat}: need at least two finished chains with mode_rb_trace.dat")
    qr = np.array([r["Q_raw"] for r in rows])
    qm = np.array([r["Q_mrb"] for r in rows])
    n = len(rows)
    ratio = qr.var(ddof=1) / qm.var(ddof=1)
    # paired: both estimators from the same chains -> variance of the difference is small; use the
    # F interval only as a conservative guide
    lo, hi = ratio / stats.f.ppf(0.95, n - 1, n - 1), ratio / stats.f.ppf(0.05, n - 1, n - 1)
    jr = np.array([r["se2_ratio"] for r in rows])
    cost = np.nanmean([r["user_s"] for r in rows]) / np.nanmean([r["user_plain_s"] for r in rows])
    out = dict(
        material=mat,
        chains=n,
        per_chain=rows,
        pooled=dict(
            Q_raw=float(qr.mean()),
            se_Q_raw=float(qr.std(ddof=1) / math.sqrt(n)),
            Q_mrb=float(qm.mean()),
            se_Q_mrb=float(qm.std(ddof=1) / math.sqrt(n)),
            max_abs_Q_diff=float(np.max(np.abs(qm - qr))),
            between_chain_var_ratio_raw_over_mrb=float(ratio),
            between_chain_var_ratio_90CI_unpaired=[float(lo), float(hi)],
            paired_jackknife_se2_ratio_mean=float(jr.mean()),
            paired_jackknife_se2_ratio_se=float(jr.std(ddof=1) / math.sqrt(n)),
            cost_ratio_user=float(cost),
            net_efficiency_jackknife=float(jr.mean() / cost),
            net_efficiency_between=float(ratio / cost),
            all_selfchecks_clean=all(r["mode_rb_mismatch"] == "0" for r in rows),
        ),
    )
    return out


def main(argv: list[str] | None = None) -> int:
    ap = make_parser(__name__, __doc__)
    ap.add_argument("runs_mrb", help="runs made with --mode-rb (mode_rb_trace.dat)")
    ap.add_argument("runs_plain", help="hook-free native runs with the same seeds (cost baseline)")
    ap.add_argument("material")
    ap.add_argument("e_bare", type=float, help="bare band energy in eV")
    ap.add_argument("tau_max", type=float, nargs="?", default=TAU_MAX_50K, help="1/kT in 1/eV (default: 50 K)")
    a = ap.parse_args(argv)
    return emit(lambda: analyze(a.runs_mrb, a.runs_plain, a.material, a.e_bare, a.tau_max))


if __name__ == "__main__":
    sys.exit(main())
