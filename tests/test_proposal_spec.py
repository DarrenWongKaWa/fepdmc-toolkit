"""M2: ProposalSpec — symmetric vs provided log-ratio; reverse is negation."""

from __future__ import annotations

import math

import pytest

from fepdmc_toolkit.diagram_compiler.proposal import ProposalError, ProposalSpec


def test_symmetric_requires_zero_log_ratio():
    spec = ProposalSpec.symmetric()
    assert spec.log_q_reverse_minus_forward(0.0) == 0.0
    with pytest.raises(ProposalError):
        spec.log_q_reverse_minus_forward(0.1)


def test_provided_ratio_is_reciprocal_under_sign_flip():
    spec = ProposalSpec.provided_log_ratio()
    value = math.log(0.3) - math.log(0.7)
    assert spec.log_q_reverse_minus_forward(value) == value
    assert spec.reverse_log_ratio(value) == pytest.approx(-value)


def test_nonfinite_log_ratio_fails_closed():
    spec = ProposalSpec.provided_log_ratio()
    with pytest.raises(ProposalError):
        spec.log_q_reverse_minus_forward(float("nan"))
    with pytest.raises(ProposalError):
        spec.log_q_reverse_minus_forward(float("inf"))


def test_proposal_spec_json_round_trip():
    for spec in (ProposalSpec.symmetric(), ProposalSpec.provided_log_ratio()):
        assert ProposalSpec.from_json(spec.to_json()) == spec
