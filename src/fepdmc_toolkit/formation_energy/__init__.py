"""LiF-electron formation energy postprocessing.

Q = Re(sum_b N_b / tauMax_b / sum_b D_b) - E_bare

The unshifted ratio (~9.105 eV on the frozen Method0 setting) is NOT Q.
"""

from __future__ import annotations

import numpy as np

E_BARE_LIF_METHOD0 = 9.35487318746596053
# 1% of |Q| ~ 0.25 eV is ~2.5 meV, not ~91 meV.
ONE_PCT_OF_QUARTER_EV_MEV = 2.5


def formation_energy_eV(ratio_eV: float, e_bare: float = E_BARE_LIF_METHOD0) -> float:
    return float(ratio_eV) - float(e_bare)


def relative_uncertainty(q_eV: float, se_eV: float) -> float:
    q = float(q_eV)
    se = float(se_eV)
    if not np.isfinite(q) or not np.isfinite(se) or se < 0.0:
        raise ValueError("Q must be finite and SE finite and >= 0")
    denom = abs(q)
    if denom == 0.0:
        raise ValueError("relative uncertainty undefined for Q=0")
    return se / denom
