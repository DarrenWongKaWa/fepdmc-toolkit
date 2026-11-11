"""C5 Q must be formation energy, not the unsubtracted ratio."""

import csv
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
FROZEN = ROOT / "benchmarks/c5_dev/frozen_results.json"


def test_reported_q_is_negative_formation_energy():
    rec = json.loads(FROZEN.read_text())
    assert "E_bare" in rec["definition"]
    bare = rec["E_bare_eV"]
    arms = rec["arms"]
    assert abs(bare - 9.35487318746596053) < 1e-12
    for arm in ("B0", "Bbest", "P1"):
        q = arms[arm]["Q_eV"]
        assert q < 0
        assert abs(q + 0.25) < 0.01


def test_hac_failures_are_recorded():
    csv_path = ROOT / "benchmarks/c5_dev/c5_per_chain_diagnostics.csv"
    rec = json.loads(FROZEN.read_text())
    assert rec["n_hac_fail"] == 5
    assert rec["hac_fail_chains"] == ["B0_r3", "Bbest_r1", "Bbest_r2", "Bbest_r4", "P1_r1"]
    rows = list(csv.DictReader(csv_path.open()))
    failed = [r["chain"] for r in rows if r["hac_pass_1.25"] in ("False", "false")]
    assert failed == ["B0_r3", "Bbest_r1", "Bbest_r2", "Bbest_r4", "P1_r1"]
