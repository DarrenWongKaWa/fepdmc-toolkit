"""M4/M5/M8: compile delayed-acceptance kernel from exact+cheap evaluators."""

from __future__ import annotations

import math

import pytest

from fepdmc_toolkit.diagram_compiler import (
    Binding,
    build_scalar_ir,
    build_twoband_ir,
    compile_cheap,
    compile_delayed_acceptance,
    compile_exact,
)
from fepdmc_toolkit.diagram_compiler.da_kernel import DAError
from fepdmc_toolkit.diagram_compiler.evaluator import torus_point
from fepdmc_toolkit.diagram_compiler.proposal import ProposalSpec


def _binding(n: int, *, L: int = 4, variant: str = "x") -> Binding:
    if variant == "x":
        k_idx, q0, tau_step, t, omega, g = 0, 1, 0.05, 1.0, 0.8, 0.5
    else:
        k_idx, q0, tau_step, t, omega, g = 1, 2, 0.07, 1.1, 0.9, 0.4
    return Binding(
        k=torus_point(k_idx, L),
        q=tuple(torus_point(q0 + i, L) for i in range(n)),
        tau=tuple(tau_step * (i + 1) for i in range(2 * n)),
        t=t,
        omega=omega,
        g=g,
        L=L,
        delta=0.35,
        gap=0.5,
    )


def test_kernel_spec_round_trips():
    ir = build_twoband_ir(1)
    kernel = compile_delayed_acceptance(ir)
    revived = type(kernel.spec).from_json(kernel.spec.to_json())
    assert revived == kernel.spec
    assert revived.design == "B"
    assert revived.cheap_policy_name == "propagator_only_v1"
    assert revived.target_policy_name == "positive_real_F_v1"


def test_twoband_exact_differs_from_cheap_and_both_positive():
    ir = build_twoband_ir(1)
    x = _binding(1)
    exact = compile_exact(ir).evaluate(x)["F"]
    cheap = compile_cheap(ir).evaluate(x)["F_hat"]
    assert exact.real > 0 and cheap.real > 0
    assert abs(exact - cheap) > 1e-8


def test_stage1_reject_skips_exact_y(monkeypatch):
    ir = build_twoband_ir(1)
    kernel = compile_delayed_acceptance(ir)
    x, y = _binding(1), _binding(1, variant="y")
    calls = {"exact": []}
    original = kernel._exact.evaluate

    def spy(binding):
        calls["exact"].append(binding)
        return original(binding)

    monkeypatch.setattr(kernel._exact, "evaluate", spy)
    ell_hat = kernel.cheap_ell_hat(x, y)
    u1 = 0.999 if ell_hat < 0 else math.exp(min(0.0, ell_hat)) + 0.1
    u1 = min(u1, 0.999999)
    # If ell_hat >= 0, stage1 always passes. Pick the direction with ell_hat < 0.
    if ell_hat >= 0:
        x, y = y, x
        ell_hat = kernel.cheap_ell_hat(x, y)
        assert ell_hat < 0
        u1 = 0.999
    result = kernel.evaluate_transition(x, y, 0.0, u1, 0.1)
    assert result.stage1_pass is False
    assert result.accepted is False
    assert result.exact_y_evaluated is False
    assert all(b != y for b in calls["exact"])


def test_stage2_can_reject_after_stage1_pass():
    ir = build_twoband_ir(1)
    kernel = compile_delayed_acceptance(ir)
    x, y = _binding(1), _binding(1, variant="y")
    probe = kernel.evaluate_transition(x, y, 0.0, 1e-16, 1e-16)
    assert probe.stage1_pass is True
    delta = probe.ell_R - probe.ell_hat
    alpha2 = 1.0 if delta >= 0.0 else math.exp(delta)
    assert alpha2 < 1.0
    u2 = (alpha2 + 1.0) / 2.0
    result = kernel.evaluate_transition(x, y, 0.0, 1e-16, u2)
    assert result.stage1_pass is True
    assert result.exact_y_evaluated is True
    assert result.stage2_pass is False
    assert result.accepted is False


def test_stage1_pass_evaluates_exact_y(monkeypatch):
    ir = build_twoband_ir(1)
    kernel = compile_delayed_acceptance(ir)
    x, y = _binding(1), _binding(1, variant="y")
    calls = []
    original = kernel._exact.evaluate

    def spy(binding):
        calls.append(binding)
        return original(binding)

    monkeypatch.setattr(kernel._exact, "evaluate", spy)
    result = kernel.evaluate_transition(x, y, 0.0, 1e-16, 1e-16)
    assert result.stage1_pass is True
    assert result.exact_y_evaluated is True
    assert any(b == y for b in calls)


def test_cached_exact_x_skips_reevaluation(monkeypatch):
    ir = build_twoband_ir(1)
    kernel = compile_delayed_acceptance(ir)
    x, y = _binding(1), _binding(1, variant="y")
    log_px = kernel.exact_log_weight(x)
    calls = []
    original = kernel._exact.evaluate

    def spy(binding):
        calls.append(binding)
        return original(binding)

    monkeypatch.setattr(kernel._exact, "evaluate", spy)
    result = kernel.evaluate_transition(
        x,
        y,
        0.0,
        1e-16,
        1e-16,
        exact_log_weight_x=log_px,
        exact_x=x,
    )
    assert result.stage1_pass is True
    assert result.exact_x_evaluated is False
    assert any(b == y for b in calls)
    assert all(b != x for b in calls)


def test_stale_exact_x_cache_fails_closed():
    ir = build_twoband_ir(1)
    kernel = compile_delayed_acceptance(ir)
    x, y = _binding(1), _binding(1, variant="y")
    with pytest.raises(DAError, match="cache"):
        kernel.evaluate_transition(
            x,
            y,
            0.0,
            1e-16,
            1e-16,
            exact_log_weight_x=0.0,
            exact_x=y,
        )


def test_g_zero_target_fails_closed():
    ir = build_scalar_ir(1)
    kernel = compile_delayed_acceptance(ir)
    x = _binding(1, variant="x")
    y = Binding(**{**x.__dict__, "g": 0.0})
    with pytest.raises((DAError, ValueError)):
        kernel.evaluate_transition(x, y, 0.0, 1e-16, 1e-16)


def test_uniform_outside_unit_interval_fails_closed():
    ir = build_scalar_ir(1)
    kernel = compile_delayed_acceptance(ir)
    x, y = _binding(1), _binding(1, variant="y")
    with pytest.raises(DAError):
        kernel.evaluate_transition(x, y, 0.0, -0.1, 0.5)
    with pytest.raises(DAError):
        kernel.evaluate_transition(x, y, 0.0, 0.5, 1.1)


def test_symmetric_proposal_rejects_nonzero_log_q():
    ir = build_scalar_ir(1)
    kernel = compile_delayed_acceptance(ir, proposal_spec=ProposalSpec.symmetric())
    x, y = _binding(1), _binding(1, variant="y")
    with pytest.raises(DAError):
        kernel.evaluate_transition(x, y, 0.3, 0.1, 0.1)
