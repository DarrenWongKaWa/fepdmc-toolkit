"""Exact-corrected delayed-acceptance algebra for analytic P1.

Native swap target is abs(Re M), not abs(complex M).
Stage-1 uses a clipped cheap score; stage-2 restores the exact native ratio.
"""

from __future__ import annotations

import math

import numpy as np

LOG10 = math.log(10.0)


def native_re_target(mat, external=1.0) -> float:
    """Native weight used in the swap ratio: Re(external * mat)."""
    z = np.asarray(external, dtype=np.complex128) * np.asarray(mat, dtype=np.complex128)
    return float(np.real(z))


def abs_complex_target(mat, external=1.0) -> float:
    """Wrong control target: |external * mat|."""
    z = np.asarray(external, dtype=np.complex128) * np.asarray(mat, dtype=np.complex128)
    return float(np.abs(z))


def native_re_ratio(mat_new, mat_old, external_new=1.0, external_old=1.0) -> float:
    """R = abs(Re M_new) / abs(Re M_old). Zero/nonfinite old is rejected."""
    py = abs(native_re_target(mat_new, external_new))
    px = abs(native_re_target(mat_old, external_old))
    if py == 0.0 and np.isfinite(px) and px > 0.0:
        return 0.0
    if px == 0.0 or not np.isfinite(px) or not np.isfinite(py):
        return float("nan")
    return py / px


def clip_log(ell: float, bound: float = LOG10) -> float:
    return float(np.clip(ell, -bound, bound))


def delayed_acceptance_prob(pi_old: float, pi_new: float, rhat: float) -> float:
    """Two-stage DA with exact correction.

    a = min(1, rhat) * min(1, R / rhat) with R = pi_new / pi_old
    and pi = abs(Re M) (passed in already as non-negative reals).
    Negative weights are errors, not negative probabilities.
    """
    px, py, rh = float(pi_old), float(pi_new), float(rhat)
    if not np.isfinite(px) or not np.isfinite(py) or not np.isfinite(rh):
        raise ValueError("DA weights and rhat must be finite")
    if px < 0.0 or py < 0.0:
        raise ValueError("DA weights must be >= 0")
    if py == 0.0 and px > 0.0:
        return 0.0
    if px <= 0.0:
        raise ValueError("DA old weight must be > 0 when new weight is positive")
    if rh <= 0.0:
        raise ValueError("rhat must be > 0")
    r = py / px
    return min(1.0, rh) * min(1.0, r / rh)


def stage2_exact_accept(ran: float, ell_r: float, ell_hat: float) -> bool:
    """Native second-stage accept: log(ran) < min(0, ellR - ell_hat).

    Matches update_swap: log(max(ran, 1e-300)) < min(0, ellR - ell_hat).
    Dropping `- ell_hat` is a different algorithm.
    """
    u, er, eh = float(ran), float(ell_r), float(ell_hat)
    if not np.isfinite(u) or u <= 0.0:
        raise ValueError("stage-2 RNG draw must be finite and > 0")
    if not np.isfinite(er) or not np.isfinite(eh):
        raise ValueError("ellR and ell_hat must be finite")
    return math.log(max(u, 1.0e-300)) < min(0.0, er - eh)


def occupancy_reverse_row(row: dict) -> dict:
    """Reverse a swap at fixed time order: swap phonon occupancy, not Δt and ΔE.

    Do not flip both (tauR-tauL) and (ek_new-ek_old); that double-flips the
    electron-propagator exponent and fails antisymmetry.
    """
    rev = dict(row)
    rev["w1"], rev["w2"] = row["w2"], row["w1"]
    rev["tpl"], rev["tpr"] = row["tpr"], row["tpl"]
    return rev


def logp_pkchange(row: dict, ek_old: float, ek_new: float, clip: float = LOG10) -> float:
    """Cheap P1 log-score from times and frequencies (no g / environment)."""

    def ld(w, d):
        w = float(w)
        if (not math.isfinite(w)) or w < 1.0e-3:
            return math.log(1.0e-15)
        return -w * abs(float(d))

    tau12 = float(row["tauR"]) - float(row["tauL"])
    ell = -(float(ek_new) - float(ek_old)) * tau12
    ell += ld(row["w1"], float(row["tauR"]) - float(row["tpl"])) - ld(row["w1"], float(row["tauL"]) - float(row["tpl"]))
    ell += ld(row["w2"], float(row["tauL"]) - float(row["tpr"])) - ld(row["w2"], float(row["tauR"]) - float(row["tpr"]))
    return clip_log(ell, clip)


def reverse_logp(row: dict, ek_old: float, ek_new: float) -> float:
    return logp_pkchange(occupancy_reverse_row(row), ek_new, ek_old)
