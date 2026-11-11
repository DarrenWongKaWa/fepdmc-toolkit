"""Formation energy is the E_bare-shifted observable, not the ~9.105 eV ratio."""

from __future__ import annotations

import json
from pathlib import Path

from fepdmc_toolkit.formation_energy import E_BARE_LIF_METHOD0, formation_energy_eV

ROOT = Path(__file__).resolve().parents[1]
FROZEN = ROOT / "benchmarks/c5_dev/frozen_results.json"


def test_e_bare_constant():
    assert abs(E_BARE_LIF_METHOD0 - 9.35487318746596053) < 1e-15


def test_unshifted_ratio_is_not_q():
    ratio = 9.104619927424677
    q = formation_energy_eV(ratio)
    assert q < 0.0
    assert abs(q + 0.25025326004128345) < 1e-12
    assert abs(ratio - 9.105) < 0.01
    assert abs(q - 9.105) > 9.0


def test_c5_frozen_q_is_negative_formation_energy():
    rec = json.loads(FROZEN.read_text())
    assert rec["E_bare_eV"] == E_BARE_LIF_METHOD0 or abs(rec["E_bare_eV"] - E_BARE_LIF_METHOD0) < 1e-12
    for vals in rec["arms"].values():
        assert vals["Q_eV"] < 0.0
        assert abs(vals["Q_eV"] + 0.25) < 0.01
        assert abs(vals["Q_eV"] - formation_energy_eV(vals["ratio_eV"])) < 1e-12
