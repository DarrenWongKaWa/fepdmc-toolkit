"""R1 scoped oracles: reuse is real; extra CSE gain was not established."""

from __future__ import annotations

import json
from pathlib import Path

from fepdmc_toolkit.r1_fixed_order import OBJECTS, scoped_status

ROOT = Path(__file__).resolve().parents[1]
FROZEN = json.loads((ROOT / "benchmarks/r1_fixed_order/frozen_results.json").read_text())


def test_three_objects_remain_distinct():
    assert OBJECTS == ("BOUND_C", "SHARED_X_GROUP", "SUMMED_KERNEL")


def test_fixed_order_r1_cse_near_unity():
    assert FROZEN["r1_over_cse_overall"] == 0.9992
    lo, hi = FROZEN["r1_over_cse_tier_range"]
    assert lo < 1.0 < hi or abs(hi - 1.0) < 0.01
    assert FROZEN["status"].startswith("VALIDATED_IN_SCOPE")
    assert "NO_ADDITIONAL_CSE_GAIN" in FROZEN["status"]


def test_material_consumer_is_future_work_not_a_universal_negative():
    rec = scoped_status()
    assert rec["material_grouped_consumer"] == "FUTURE_WORK"
    assert rec["additional_cse_gain"] == "NO_ADDITIONAL_CSE_GAIN_ESTABLISHED"
    assert FROZEN["material_consumer"] == "FUTURE_WORK"
    assert rec["evaluator_in_this_package"] is False
    assert rec["public_release"] == "HISTORICAL_RESULT_METADATA"


def test_public_module_has_no_evaluator_entry_points():
    import fepdmc_toolkit.r1_fixed_order as mod

    banned = ("evaluate", "compile", "rebind", "build_plan", "cse")
    names = set(dir(mod))
    assert not any(b in names for b in banned)


def test_abs_sum_is_not_sum_abs():
    # Finite-domain fact that blocks naive grouped |sum D| sampling.
    d = [1.0, -1.0, 0.5]
    assert abs(sum(d)) != sum(abs(x) for x in d)
