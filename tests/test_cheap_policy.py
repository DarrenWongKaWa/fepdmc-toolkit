"""M1: CheapPolicy schema — KEEP / DROP / APPROXIMATE / REFUSE."""

from __future__ import annotations

import json

import pytest

from fepdmc_toolkit.diagram_compiler.cheap_policy import (
    POLICY_PROPAGATOR_ONLY_V1,
    CheapPolicy,
    CheapPolicyError,
    propagator_only_v1,
)


def test_propagator_only_v1_lists_are_disjoint_and_versioned():
    policy = propagator_only_v1()
    assert policy.schema_version == 1
    assert policy.name == POLICY_PROPAGATOR_ONLY_V1
    listed = policy.keep + policy.drop + policy.approximate + policy.refuse
    assert len(listed) == len(set(listed))
    assert "factor.electron_propagator" in policy.keep
    assert "factor.phonon_propagator" in policy.keep
    assert "factor.vertex" in policy.drop
    assert "electron_interface.matrix2" in policy.approximate
    assert "eval.matmul" in policy.drop
    assert "eval.trace" in policy.drop


def test_policy_json_round_trip_is_deterministic():
    policy = propagator_only_v1()
    text = policy.to_json()
    revived = CheapPolicy.from_json(text)
    assert revived == policy
    again = json.loads(revived.to_json())
    assert list(again.keys()) == sorted(again.keys())
    assert again == json.loads(text)


def test_unknown_category_is_uncovered_not_silently_kept():
    policy = propagator_only_v1()
    with pytest.raises(CheapPolicyError, match="uncovered"):
        policy.action("factor.dressed_propagator")


def test_refuse_category_is_explicit():
    policy = CheapPolicy(
        schema_version=1,
        name="hand_refuse",
        keep=("factor.electron_propagator",),
        drop=(),
        approximate=(),
        refuse=("factor.vertex",),
    )
    assert policy.action("factor.vertex") == "refuse"
    with pytest.raises(CheapPolicyError, match="refuse"):
        policy.require_allowed("factor.vertex")


def test_duplicate_category_in_two_lists_is_rejected():
    with pytest.raises(CheapPolicyError, match="duplicate"):
        CheapPolicy(
            schema_version=1,
            name="bad",
            keep=("factor.vertex",),
            drop=("factor.vertex",),
            approximate=(),
            refuse=(),
        )


def test_checked_in_policy_artifact_round_trips():
    from pathlib import Path

    root = Path(__file__).resolve().parents[1]
    text = (root / "validation/diagram_compiler/CHEAP_POLICY.json").read_text()
    assert CheapPolicy.from_json(text) == propagator_only_v1()


def test_from_dict_rejects_unknown_fields_and_functions():
    payload = propagator_only_v1().to_dict()
    payload["callback"] = "not-allowed"
    with pytest.raises(CheapPolicyError):
        CheapPolicy.from_dict(payload)
