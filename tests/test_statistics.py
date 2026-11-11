"""Frozen protocol: HAC, relative precision, P1 classification. No post-hoc retune."""

from __future__ import annotations

import json
from pathlib import Path

from fepdmc_toolkit.stats import classify_p1_wall_ratio, hac_pass

ROOT = Path(__file__).resolve().parents[1]
C5 = json.loads((ROOT / "benchmarks/c5_dev/frozen_results.json").read_text())
P1 = json.loads((ROOT / "benchmarks/p1_vs_bbest/frozen_results.json").read_text())


def test_c5_relative_precision_is_few_percent_not_one_percent_certification():
    p1 = C5["arms"]["P1"]
    assert 0.009 < p1["relative_JK_SE"] < 0.010
    assert 0.024 < p1["relative_point_95_halfwidth"] < 0.025
    # ~2.4% point half-width is not a 1% ground-state certificate
    assert p1["relative_point_95_halfwidth"] > 0.01
    assert C5["not_1pct_ground_state"] is True


def test_c5_hac_failures_not_rewritten():
    assert C5["hac_fail_chains"] == [
        "B0_r3",
        "Bbest_r1",
        "Bbest_r2",
        "Bbest_r4",
        "P1_r1",
    ]
    assert C5["n_hac_fail"] == 5


def test_predeclared_hac_threshold_unchanged():
    assert hac_pass([1.0, 1.24, 1.24]) is True
    assert hac_pass([1.0, 1.26, 1.0]) is False


def test_clean_p1_classification_crosses_threshold():
    rec = P1
    assert rec["T_boot_2.5"] < 0.95 < rec["T_boot_97.5"]
    assert rec["hac_ok"] is False
    assert classify_p1_wall_ratio(rec["T_boot_2.5"], rec["T_boot_97.5"], rec["hac_ok"]) == "P1_UNRESOLVED_WITHIN_BUDGET"
    assert rec["classification"] == "P1_UNRESOLVED_WITHIN_BUDGET"
    # K is reported; it is not a second frozen pass/fail.
    assert rec["K_ratio_point"] > 1.0
    assert rec["K_boot_2.5"] < 1.0 < rec["K_boot_97.5"]


def test_counter_identity_on_confirm_sums():
    tot = P1["p1_counters_sum"]
    assert tot["n_stage1_rej"] + tot["n_stage2"] + tot["n_accept"] == tot["n_eligible"]
    assert tot["n_exact_stage"] == tot["n_eligible"] - tot["n_stage1_rej"]
    assert 0.33 < tot["a1_rate"] < 0.34
