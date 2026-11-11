"""Correlated-MC diagnostics used in the frozen P1 protocol.

HAC Bartlett lags 2/4/8 with maxSE/minSE <= 1.25 was predeclared.
Do not retune lags after seeing results.

`mean_jackknife_se` is an arithmetic-mean jackknife. It is not the
project's ratio-of-sums jackknife. Use `ratio_of_sums_jackknife_se`
for N/D pooled ratios.
"""

from __future__ import annotations

from collections.abc import Iterable, Sequence

import numpy as np

T_CRIT_DF5_975 = 2.5705818366147395  # Student-t, two-sided 95%, df=5 (6 chains)
HAC_STABILITY_MAX = 1.25
PRACTICAL_GAIN_THRESHOLD = 0.95


def _finite_nonneg_se(se: float, name: str = "se") -> float:
    x = float(se)
    if not np.isfinite(x) or x < 0.0:
        raise ValueError(f"{name} must be finite and >= 0, got {se!r}")
    return x


def mean_jackknife_se(values: Sequence[float]) -> float:
    """Delete-one jackknife SE of an arithmetic mean. Not a ratio JK."""
    x = np.asarray(values, dtype=np.float64)
    if x.size < 2 or not np.all(np.isfinite(x)):
        raise ValueError("mean_jackknife_se requires >=2 finite values")
    n = x.size
    mean = x.mean()
    leave = np.array([(n * mean - x[i]) / (n - 1) for i in range(n)])
    return float(np.sqrt((n - 1) / n * np.sum((leave - leave.mean()) ** 2)))


def jackknife_se(values: Sequence[float]) -> float:
    """Deprecated alias of mean_jackknife_se. Not ratio-of-sums JK."""
    return mean_jackknife_se(values)


def ratio_of_sums_jackknife_se(numerators: Sequence[float], denominators: Sequence[float]) -> float:
    """Delete-one jackknife SE of (sum N)/(sum D)."""
    n = np.asarray(numerators, dtype=np.float64)
    d = np.asarray(denominators, dtype=np.float64)
    if n.shape != d.shape or n.size < 2:
        raise ValueError("N and D must be same length >= 2")
    if not np.all(np.isfinite(n)) or not np.all(np.isfinite(d)):
        raise ValueError("N and D must be finite")
    if np.any(d == 0.0) or (d.sum() - d == 0.0).any():
        raise ValueError("D and leave-one-out D sums must be nonzero")
    ntot, dtot = n.sum(), d.sum()
    leave = (ntot - n) / (dtot - d)
    k = n.size
    return float(np.sqrt((k - 1) / k * np.sum((leave - leave.mean()) ** 2)))


def point_halfwidth_t975(se: float, t_crit: float = T_CRIT_DF5_975) -> float:
    return _finite_nonneg_se(se) * float(t_crit)


def hac_pass(se_lags: Iterable[float], max_ratio: float = HAC_STABILITY_MAX) -> bool:
    """True only if every provided window is finite and positive and stable.

    Invalid windows are not dropped. NaN/Inf/zero/negative raise.
    """
    se = list(se_lags)
    if len(se) < 2:
        raise ValueError("HAC requires at least two lag windows")
    cleaned = []
    for s in se:
        x = float(s)
        if not np.isfinite(x) or x <= 0.0:
            raise ValueError(f"HAC window must be finite and > 0, got {s!r}")
        cleaned.append(x)
    return (max(cleaned) / min(cleaned)) <= float(max_ratio)


def classify_p1_wall_ratio(
    boot_low: float,
    boot_high: float,
    hac_ok: bool,
    threshold: float = PRACTICAL_GAIN_THRESHOLD,
) -> str:
    """Frozen classification. Do not edit after seeing confirmation data."""
    lo, hi, thr = float(boot_low), float(boot_high), float(threshold)
    if not np.isfinite(lo) or not np.isfinite(hi) or not np.isfinite(thr):
        raise ValueError("bootstrap endpoints and threshold must be finite")
    if lo > hi:
        raise ValueError(f"bootstrap interval reversed: {lo} > {hi}")
    if not isinstance(hac_ok, (bool, np.bool_)):
        raise ValueError("hac_ok must be bool")
    if hi < thr and bool(hac_ok):
        return "P1_GAIN_ESTABLISHED"
    if lo > thr:
        return "P1_NO_PRACTICAL_GAIN"
    return "P1_UNRESOLVED_WITHIN_BUDGET"
