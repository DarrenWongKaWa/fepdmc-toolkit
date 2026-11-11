"""Parse LINEAR_DA_COUNTS and fixture summaries.

Production counters are independent of P1_FIXTURE / dump switches.
Format: LINEAR_DA_COUNTS <mode> <score> n_invoked n_eligible n_score
        n_stage1_rej n_stage2 n_accept aux_seed
"""

from __future__ import annotations

from typing import Any


def parse_linear_da_counts(line: str) -> dict[str, Any]:
    parts = line.strip().split()
    if len(parts) != 10 or parts[0] != "LINEAR_DA_COUNTS":
        raise ValueError(f"expected 10 LINEAR_DA_COUNTS fields, got {line!r}")
    try:
        nums = [int(x) for x in parts[3:10]]
    except ValueError as exc:
        raise ValueError(f"non-integer LINEAR_DA_COUNTS fields: {line!r}") from exc
    invoked, eligible, score, a1, a2, nacc, seed = nums
    for name, val in (
        ("n_invoked", invoked),
        ("n_eligible", eligible),
        ("n_score", score),
        ("n_stage1_rej", a1),
        ("n_stage2", a2),
        ("n_accept", nacc),
    ):
        if val < 0:
            raise ValueError(f"{name} must be >= 0, got {val}")
    if invoked < eligible or eligible < score:
        raise ValueError("counter hierarchy n_invoked >= n_eligible >= n_score failed")
    if a1 + a2 + nacc != eligible:
        raise ValueError(f"counter identity failed: {line}")
    exact_stage = eligible - a1
    if exact_stage != a2 + nacc:
        raise ValueError("exact-stage identity failed")
    return {
        "mode": parts[1],
        "score": parts[2],
        "n_invoked": invoked,
        "n_eligible": eligible,
        "n_score": score,
        "n_stage1_rej": a1,
        "n_stage2": a2,
        "n_accept": nacc,
        "n_exact_stage": exact_stage,
        "aux_seed": seed,
        "a1_rate": (a1 / eligible) if eligible else None,
    }
