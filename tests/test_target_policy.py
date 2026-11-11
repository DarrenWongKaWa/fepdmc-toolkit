"""M1: TargetPolicy — positive exact F only; no silent |F| or F^2."""

from __future__ import annotations

import pytest

from fepdmc_toolkit.diagram_compiler.target_policy import (
    TargetPolicyError,
    positive_real_F_v1,
)


def test_policy_is_versioned_and_refuses_absolute_value():
    policy = positive_real_F_v1()
    assert policy.schema_version == 1
    assert policy.name == "positive_real_F_v1"
    assert policy.target == "exact_F"
    assert "abs_F" not in policy.to_dict().values()
    assert policy.to_dict()["refuse_transforms"] == ["abs_F", "F_squared", "complex_modulus"]


def test_positive_real_F_is_accepted():
    policy = positive_real_F_v1()
    log_w = policy.log_weight(0.25 + 0j)
    assert log_w == pytest.approx(__import__("math").log(0.25))


def test_zero_negative_complex_nan_inf_fail_closed():
    policy = positive_real_F_v1()
    with pytest.raises(TargetPolicyError, match="non-positive"):
        policy.log_weight(0.0 + 0j)
    with pytest.raises(TargetPolicyError, match="non-positive"):
        policy.log_weight(-0.1 + 0j)
    with pytest.raises(TargetPolicyError, match="complex"):
        policy.log_weight(0.2 + 0.01j)
    with pytest.raises(TargetPolicyError):
        policy.log_weight(float("nan") + 0j)
    with pytest.raises(TargetPolicyError):
        policy.log_weight(float("inf") + 0j)


def test_target_policy_json_round_trip():
    policy = positive_real_F_v1()
    revived = type(policy).from_json(policy.to_json())
    assert revived == policy
    extra = policy.to_dict()
    extra["callback"] = "nope"
    with pytest.raises(TargetPolicyError):
        type(policy).from_dict(extra)
