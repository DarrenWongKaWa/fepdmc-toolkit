"""Recorded result of the fixed-order R1 evaluator benchmark (metadata only).

This module is not an R1 evaluator. The recursive R1 builder used as an oracle
by the diagram compiler is in ``fepdmc_toolkit.scheme_d.physics_first.symbolic``;
the timing comparison against an optimized common-subexpression baseline is
recorded here and in benchmarks/r1_fixed_order/.
"""

from __future__ import annotations

R1_OVER_CSE_OVERALL = 0.9992
R1_OVER_CSE_TIER_MIN = 0.9863
R1_OVER_CSE_TIER_MAX = 1.0034

OBJECTS = ("BOUND_C", "SHARED_X_GROUP", "SUMMED_KERNEL")
PUBLIC_RELEASE = "HISTORICAL_RESULT_METADATA"


def scoped_status() -> dict:
    return {
        "public_release": PUBLIC_RELEASE,
        "evaluator_in_this_package": False,
        "evaluator": "VALIDATED_IN_SCOPE",
        "additional_cse_gain": "NO_ADDITIONAL_CSE_GAIN_ESTABLISHED",
        "r1_over_cse_overall": R1_OVER_CSE_OVERALL,
        "material_grouped_consumer": "FUTURE_WORK",
        "objects": OBJECTS,
    }
