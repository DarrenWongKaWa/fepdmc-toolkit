"""How much of the native Q variance is the sign? (per-measurement, blocking tau).

h   = (s O - R s) / <s>        the ratio-estimator influence of the native chain
h_O = (O - <O>)                the same chain with the sign set to +1
If sigma^2(h_O) ~ sigma^2(h), removing the sign cannot help Q much, whatever
1/<s>^2 says: the sign doubles the static variance but, being fast noise, it
also decorrelates the slow O fluctuations.
"""

import glob
import sys

import numpy as np

from ._common import TAU_MAX_50K, AnalysisError, emit, make_parser


def tint(x):
    x = np.asarray(x, float)
    v0 = x.var()
    best = 1.0
    b = 1
    while len(x) // b >= 30:
        n = len(x) // b
        m = x[: n * b].reshape(n, b).mean(1)
        best = max(best, b * m.var() / v0)
        b *= 2
    return float(best)


def analyze(runs: str, mat: str, tau: float = TAU_MAX_50K) -> dict:
    rows = []
    for f in sorted(glob.glob(f"{runs}/{mat}/chain*/rb_trace.dat")):
        a = np.loadtxt(f)
        num, den = a[:, 0], a[:, 1]
        O = num / den
        R = num.sum() / den.sum()
        s = den.mean()
        h = (num - R * den) / s / tau
        hO = (O - O.mean()) / tau
        rows.append(
            dict(
                chain=f.split("/")[-2],
                sign=float(s),
                static_h=float(h.var()),
                tau_h=tint(h),
                sigma2_h=float(h.var() * tint(h)),
                static_O=float(hO.var()),
                tau_O=tint(hO),
                sigma2_O=float(hO.var() * tint(hO)),
                tau_order=tint(a[:, 6]),
            )
        )
    if not rows:
        raise AnalysisError(f"{runs}/{mat}: no chain with rb_trace.dat")
    keys = [k for k in rows[0] if k != "chain"]
    pooled = {k: float(np.mean([r[k] for r in rows])) for k in keys}
    pooled["fraction_of_variance_not_from_sign"] = pooled["sigma2_O"] / pooled["sigma2_h"]
    pooled["Q_variance_ceiling_from_removing_sign"] = pooled["sigma2_h"] / pooled["sigma2_O"]
    pooled["naive_sign_ceiling_1_over_s2"] = 1 / pooled["sign"] ** 2
    return dict(material=mat, per_chain=rows, pooled=pooled)


def main(argv: list[str] | None = None) -> int:
    ap = make_parser(__name__, __doc__)
    ap.add_argument("runs", help="runs with rb_trace.dat")
    ap.add_argument("material")
    ap.add_argument("tau_max", type=float, nargs="?", default=TAU_MAX_50K, help="1/kT in 1/eV (default: 50 K)")
    a = ap.parse_args(argv)
    return emit(lambda: analyze(a.runs, a.material, a.tau_max))


if __name__ == "__main__":
    sys.exit(main())
