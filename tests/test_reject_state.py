"""Archived fixture *record* checks. Does not execute native rollback."""

from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def _fixture():
    for path in (ROOT / "benchmarks/p1_vs_bbest/frozen_results.json",):
        if path.exists():
            rec = json.loads(path.read_text())
            return rec.get("fixture", rec)
    raise AssertionError("fixture summary missing")


def test_stage1_rejects_recorded_zero_g_env():
    rec = _fixture()
    assert rec["ok"] is True
    assert rec["n_a1"] == 4372
    assert rec["a1_fp_match"] == 4372
    assert rec["a1_zero_g_env"] is True


def test_committed_y_p1_reverse():
    rec = _fixture()
    assert rec["n_accept_pairs"] == 6350
    assert rec["y_reverse"] == 6350


def test_dump_independent_stage1_counter_in_fixture():
    rec = _fixture()
    # LINEAR_DA_COUNTS n_invoked n_eligible n_score n_stage1_rej
    counts = rec["counts_line"]
    assert int(counts[3]) == rec["n_a1"]
