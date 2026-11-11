"""Summary-consistency replay of published JSON/CSV (not raw N/D/HAC/bootstrap).

Does not rebuild jackknife or HAC from Monte Carlo blocks. Clean P1
confirmation has no public per-block N/D series.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from fepdmc_toolkit.formation_energy import E_BARE_LIF_METHOD0, formation_energy_eV  # noqa: E402
from fepdmc_toolkit.stats import classify_p1_wall_ratio  # noqa: E402


def main() -> int:
    print("MODE=summary-consistency; not raw-stat reconstruction")
    c5 = json.loads((ROOT / "benchmarks/c5_dev/frozen_results.json").read_text())
    p1 = json.loads((ROOT / "benchmarks/p1_vs_bbest/frozen_results.json").read_text())
    print("E_bare =", E_BARE_LIF_METHOD0)
    print("C5 DEV formation energies (fixed setting, not 1% GS):")
    for name, arm in c5["arms"].items():
        q = formation_energy_eV(arm["ratio_eV"])
        assert abs(q - arm["Q_eV"]) < 1e-12
        print(
            f"  {name:6s}  Q={arm['Q_eV']:+.9f} eV  "
            f"relSE={100 * arm['relative_JK_SE']:.3f}%  "
            f"relHW={100 * arm['relative_point_95_halfwidth']:.3f}%"
        )
    t_p1 = p1["arms"]["P1"]["T_mean_s"]
    t_b = p1["arms"]["Bbest"]["T_mean_s"]
    t_ratio = t_p1 / t_b
    assert abs(t_ratio - p1["T_ratio_point"]) < 1e-12
    k_p1 = t_p1 * p1["arms"]["P1"]["VarJK"]
    k_b = t_b * p1["arms"]["Bbest"]["VarJK"]
    k_ratio = k_p1 / k_b
    assert abs(k_ratio - p1["K_ratio_point"]) < 1e-12
    cls = classify_p1_wall_ratio(p1["T_boot_2.5"], p1["T_boot_97.5"], p1["hac_ok"])
    print(f"clean P1 / B-best T ratio {t_ratio:.4f}  [{p1['T_boot_2.5']:.4f}, {p1['T_boot_97.5']:.4f}]  class={cls}")
    print(
        f"K=T*VarJK point {k_ratio:.3f}  [{p1['K_boot_2.5']:.3f}, {p1['K_boot_97.5']:.1f}]  DIAGNOSTIC_NOT_A_PASS_RULE"
    )
    assert cls == "P1_UNRESOLVED_WITHIN_BUDGET"
    print("ok")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
