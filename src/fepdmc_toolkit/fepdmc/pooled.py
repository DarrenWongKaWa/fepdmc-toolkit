"""Native chain A vs a modified chain C, by the pooled ratio estimator over chains.

All chains of one side sample the same measure, so the estimate is the pooled
ratio R = sum_i n_i / sum_i d_i (n_i = mean numerator, d_i = mean sign of chain i),
not the mean of per-chain ratios: in the low-sign regime the latter is biased by
chains that sit in low-sign regions. Per chain, n_i = Etrue_i * tau_max * <g/Z>_i
from stdout. Chain-level influence h_i = (n_i - R d_i) / mean(d) gives
Var(R) = var(h) / N; var_A(h)/var_C(h) is the per-step variance ratio (90% F
interval), CPU user time gives the cost ratio, efficiency = ratio / cost.

Optional trace diagnostics (--trace-a/--trace-c: files with numerator, sign and
order columns): spread of chain signs, drift of 10 block means of the order, and
blocking tau_int of the order.

The 90% F interval needs scipy; without it the interval is reported as null.

Usage: fepdmc-toolkit compare <runs_A> <runs_C> <material> <E_bare> <tau_max>
         [--trace-a DIR FILE NUMCOL SIGNCOL ORDERCOL] [--trace-c DIR FILE NUMCOL SIGNCOL ORDERCOL]
"""

from __future__ import annotations

import argparse
import json
import math
import re
import sys
from pathlib import Path

import numpy as np


def cpu(c: Path) -> float:
    """CPU user time in seconds from the ``time`` output in stderr.log; NaN if it is absent."""
    log = c / "stderr.log"
    t = re.search(r"user\s+(\d+)m([\d.]+)s", log.read_text()) if log.is_file() else None
    return int(t.group(1)) * 60 + float(t.group(2)) if t else float("nan")


def side(runs: str, mat: str, tau: float, e_bare: float) -> dict:
    """Pooled estimate and chain-level statistics of the runs in ``<runs>/<mat>/chain*``.

    Raises ValueError if fewer than two chains have finished (printed Etrue and <g/Z>).
    """
    rows = []
    for c in sorted((Path(runs) / mat).glob("chain*")):
        log = c / "stdout.log"
        s = log.read_text() if log.is_file() else ""
        e = re.search(r"Etrue =\s+(\S+)", s)
        g = re.search(r"<g/Z> =\s+(\S+)", s)
        if not (e and g):
            continue
        E, d = float(e.group(1)), float(g.group(1))
        rows.append(dict(chain=c.name, n=E * tau * d, d=d, Q_chain=E - e_bare, cpu_s=cpu(c)))
    if len(rows) < 2:
        raise ValueError(
            f"{Path(runs) / mat}: {len(rows)} finished chain(s) found; the pooled estimator needs at least two"
        )
    cpus = np.array([r["cpu_s"] for r in rows])
    n = np.array([r["n"] for r in rows])
    d = np.array([r["d"] for r in rows])
    k = len(rows)
    R = n.sum() / d.sum()
    h = (n - R * d) / d.mean()
    jk = np.array([(n.sum() - n[i]) / (d.sum() - d[i]) for i in range(k)])
    return dict(
        n_chains=k,
        Q_pooled=R / tau - e_bare,
        se_delta=float(h.std(ddof=1) / math.sqrt(k) / tau),
        se_jackknife=float(math.sqrt((k - 1) / k * ((jk - jk.mean()) ** 2).sum()) / tau),
        var_h=float(h.var(ddof=1)),
        sign_mean=float(d.mean()),
        sign_min=float(d.min()),
        sign_max=float(d.max()),
        sign_sd=float(d.std(ddof=1)),
        Q_mean_of_chain_ratios=float(np.mean([r["Q_chain"] for r in rows])),
        cpu_s=float("nan") if np.isnan(cpus).all() else float(np.nanmean(cpus)),
        per_chain=rows,
    )


def tau_int(x) -> float:
    x = np.asarray(x, float)
    v0 = x.var()
    best, b = 1.0, 1
    while len(x) // b >= 30:
        m = len(x) // b
        best = max(best, b * x[: m * b].reshape(m, b).mean(1).var() / v0)
        b *= 2
    return float(best)


def traces(spec, mat: str) -> dict:
    d, f, cn, cs, co = spec
    cn, cs, co = int(cn), int(cs), int(co)
    drift, tord, tsign = [], [], []
    for c in sorted((Path(d) / mat).glob("chain*")):
        a = np.loadtxt(c / f)
        o = a[:, co]
        s = a[:, cs]
        nb = len(o) // 10
        drift.append(float(o[: nb * 10].reshape(10, nb).mean(1).std(ddof=1)))
        tord.append(tau_int(o))
        tsign.append(tau_int(s))
    return dict(
        order_block_mean_sd=float(np.mean(drift)),
        tau_int_order=float(np.mean(tord)),
        tau_int_sign=float(np.mean(tsign)),
    )


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        prog="fepdmc-toolkit compare", description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("runs_a", help="run directory of side A (contains <material>/chainNN)")
    ap.add_argument("runs_c", help="run directory of side C")
    ap.add_argument("material", help="material subdirectory, for example lif_hole")
    ap.add_argument("e_bare", type=float, help="bare band energy in eV; Q = E - E_bare")
    ap.add_argument("tau", type=float, help="tau_max = 1/kT in 1/eV")
    trace = ("DIR", "FILE", "NUMCOL", "SIGNCOL", "ORDERCOL")
    ap.add_argument("--trace-a", nargs=5, metavar=trace, help="trace diagnostics for side A")
    ap.add_argument("--trace-c", nargs=5, metavar=trace, help="trace diagnostics for side C")
    a = ap.parse_args(argv)
    try:
        A = side(a.runs_a, a.material, a.tau, a.e_bare)
        C = side(a.runs_c, a.material, a.tau, a.e_bare)
    except ValueError as exc:
        print(exc, file=sys.stderr)
        return 2
    ratio = A["var_h"] / C["var_h"]
    try:
        from scipy import stats  # only needed for the F interval

        lo = ratio / stats.f.ppf(0.95, A["n_chains"] - 1, C["n_chains"] - 1)
        hi = ratio / stats.f.ppf(0.05, A["n_chains"] - 1, C["n_chains"] - 1)
    except ImportError:
        lo = hi = None
    cost = C["cpu_s"] / A["cpu_s"]
    out = dict(
        material=a.material,
        A=A,
        C=C,
        comparison=dict(
            Q_difference=C["Q_pooled"] - A["Q_pooled"],
            Q_difference_sigma=(C["Q_pooled"] - A["Q_pooled"]) / math.hypot(A["se_jackknife"], C["se_jackknife"]),
            variance_ratio_A_over_C=ratio,
            variance_ratio_90CI=[lo, hi],
            cpu_ratio_C_over_A=cost,
            efficiency=ratio / cost,
            efficiency_90CI=[lo / cost, hi / cost] if lo is not None else None,
        ),
    )
    if a.trace_a:
        out["A"]["traces"] = traces(a.trace_a, a.material)
    if a.trace_c:
        out["C"]["traces"] = traces(a.trace_c, a.material)
    json.dump(out, sys.stdout, indent=2)
    sys.stdout.write("\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
