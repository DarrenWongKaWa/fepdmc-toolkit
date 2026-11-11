"""M2/M4/M5: cheap lowering from the same DiagramIR, fail-closed uncovered nodes."""

from __future__ import annotations

import pytest

from fepdmc_toolkit.diagram_compiler import (
    Binding,
    build_scalar_ir,
    build_twoband_ir,
    compile_cheap,
    compile_evaluator,
    compile_exact,
    evaluate_diagramwise,
    numeric_close,
)
from fepdmc_toolkit.diagram_compiler.cheap_policy import CheapPolicy, CheapPolicyError
from fepdmc_toolkit.diagram_compiler.evaluator import torus_point
from fepdmc_toolkit.diagram_compiler.ir import KIND_ELECTRON, Factor, MomentumForm


def _binding(n: int, *, L: int = 4, g: float = 0.5, **extra) -> Binding:
    tau = tuple(0.05 * (i + 1) for i in range(2 * n))
    q = tuple(torus_point(i + 1, L) for i in range(n))
    fields = dict(
        k=torus_point(0, L),
        q=q,
        tau=tau,
        t=1.0,
        omega=0.8,
        g=g,
        L=L,
        delta=0.35,
        gap=0.5,
    )
    fields.update(extra)
    return Binding(**fields)


def test_hand_n1_scalar_cheap_matches_exact():
    ir = build_scalar_ir(1)
    binding = _binding(1)
    exact = compile_exact(ir).evaluate(binding)
    cheap = compile_cheap(ir, policy="propagator_only_v1").evaluate(binding)
    assert numeric_close(cheap["F_hat"], exact["F"])
    assert cheap["log_W_hat"] == pytest.approx(float(__import__("math").log(exact["F"].real)))


def test_compile_cheap_snapshots_ir():
    ir = build_scalar_ir(1)
    binding = _binding(1)
    cheap = compile_cheap(ir)
    before = cheap.evaluate(binding)
    ir.diagrams[0].sign = -1
    after = cheap.evaluate(binding)
    assert after == before


def test_scalar_n2_generated_cheap_runs_and_equals_exact():
    ir = build_scalar_ir(2)
    binding = _binding(2)
    cheap = compile_cheap(ir)
    exact_f = compile_evaluator(ir)(binding)["F"]
    assert numeric_close(cheap.evaluate(binding)["F_hat"], exact_f)


@pytest.mark.parametrize("n", [1, 2, 3, 4])
def test_scalar_propagator_only_equals_exact(n):
    ir = build_scalar_ir(n)
    binding = _binding(n)
    cheap = compile_cheap(ir).evaluate(binding)
    exact = evaluate_diagramwise(ir, binding)
    assert numeric_close(cheap["F_hat"], exact["F"])


@pytest.mark.parametrize("n", [1, 2, 3])
def test_twoband_cheap_equals_scalar_exact_not_twoband_exact(n):
    scalar = build_scalar_ir(n)
    twoband = build_twoband_ir(n)
    binding = _binding(n)
    cheap = compile_cheap(twoband).evaluate(binding)
    scalar_f = evaluate_diagramwise(scalar, binding)["F"]
    twoband_f = evaluate_diagramwise(twoband, binding)["F"]
    assert numeric_close(cheap["F_hat"], scalar_f)
    assert not numeric_close(cheap["F_hat"], twoband_f)


def test_twoband_cheap_has_no_op_output():
    cheap = compile_cheap(build_twoband_ir(1)).evaluate(_binding(1))
    assert "Op" not in cheap
    assert "F_hat" in cheap and "log_W_hat" in cheap


def test_incomplete_policy_fails_closed_on_present_electron():
    ir = build_scalar_ir(1)
    policy = CheapPolicy(
        schema_version=1,
        name="missing_electron",
        keep=(
            "factor.phonon_propagator",
            "factor.prefactor",
            "diagram.sign",
            "diagram.multiplicity",
            "family.global_factors",
            "diagram_sum.all_pairings",
            "momentum.q_terms",
        ),
        drop=("factor.vertex", "eval.matmul", "eval.trace"),
        approximate=(),
        refuse=(),
    )
    with pytest.raises(CheapPolicyError, match="uncovered"):
        compile_cheap(ir, policy=policy)


def test_unknown_factor_kind_fails_closed_not_treated_as_one():
    ir = build_scalar_ir(1)
    fake = Factor(KIND_ELECTRON, interval=(0, 1), open_chord_slots=(0,), k_form=MomentumForm(1, ((0, -1),)))
    fake.kind = "dressed_propagator"
    ir.diagrams[0].factors = (fake,) + ir.diagrams[0].factors[1:]
    with pytest.raises((CheapPolicyError, ValueError)):
        compile_cheap(ir)


def test_cheap_evaluator_json_round_trip():
    ir = build_scalar_ir(2)
    binding = _binding(2)
    cheap = compile_cheap(ir)
    revived = type(cheap).from_json(cheap.to_json())
    assert revived.policy == cheap.policy
    assert numeric_close(revived.evaluate(binding)["F_hat"], cheap.evaluate(binding)["F_hat"])


def test_twoband_refuses_when_matrix_interface_is_not_approximated():
    ir = build_twoband_ir(1)
    policy = CheapPolicy(
        schema_version=1,
        name="drop_vertex_only",
        keep=(
            "factor.electron_propagator",
            "factor.phonon_propagator",
            "factor.prefactor",
            "diagram.sign",
            "diagram.multiplicity",
            "family.global_factors",
            "diagram_sum.all_pairings",
            "momentum.q_terms",
        ),
        drop=("factor.vertex", "eval.matmul", "eval.trace"),
        approximate=(),
        refuse=(),
    )
    with pytest.raises(CheapPolicyError, match="uncovered"):
        compile_cheap(ir, policy=policy)


def test_prefactor_is_required_for_a_pure_g_transition():
    import math

    ir = build_scalar_ir(1)
    cheap = compile_cheap(ir)
    x = _binding(1, g=0.5)
    y = _binding(1, g=0.25)
    ell = cheap.transition_score(x, y)
    assert ell == pytest.approx(math.log((0.25 / 0.5) ** 2))


def test_dropping_prefactor_zeros_a_pure_g_transition():
    ir = build_scalar_ir(1)
    policy = CheapPolicy(
        schema_version=1,
        name="drop_prefactor",
        keep=(
            "factor.electron_propagator",
            "factor.phonon_propagator",
            "diagram.sign",
            "diagram.multiplicity",
            "diagram_sum.all_pairings",
            "momentum.q_terms",
        ),
        drop=("factor.prefactor", "family.global_factors", "factor.vertex", "eval.matmul", "eval.trace"),
        approximate=(),
        refuse=(),
    )
    cheap = compile_cheap(ir, policy=policy)
    x = _binding(1, g=0.5)
    y = _binding(1, g=0.25)
    assert cheap.transition_score(x, y) == pytest.approx(0.0)


def test_unimplemented_q_term_drop_fails_closed():
    ir = build_scalar_ir(2)
    policy = CheapPolicy(
        schema_version=1,
        name="drop_q_terms",
        keep=(
            "factor.electron_propagator",
            "factor.phonon_propagator",
            "factor.prefactor",
            "diagram.sign",
            "diagram.multiplicity",
            "family.global_factors",
            "diagram_sum.all_pairings",
        ),
        drop=("momentum.q_terms", "factor.vertex", "eval.matmul", "eval.trace"),
        approximate=(),
        refuse=(),
    )
    with pytest.raises(CheapPolicyError, match="not implemented"):
        compile_cheap(ir, policy=policy)


def test_g_zero_log_weight_fails_closed():
    ir = build_scalar_ir(1)
    cheap = compile_cheap(ir)
    with pytest.raises(ValueError, match="non-positive"):
        cheap.log_weight(_binding(1, g=0.0))
