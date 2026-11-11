"""M3/M4: differential tests of the generated evaluator against the R1
oracle path (r1_symbolic_build + ReferenceBinder + exec) and the independent
index_oracle enumerators. Tolerance is the repository numeric standard
(atol 1e-12 + rtol 1e-10), fixed before any results were examined.

Mutation tests then verify the differential harness is load-bearing: each
physics fault must move the generated value away from the R1 oracle.
"""

from __future__ import annotations

import numpy as np
import pytest

from fepdmc_toolkit.diagram_compiler import (
    Diagram,
    DiagramIR,
    Factor,
    MomentumForm,
    evaluate_diagramwise,
    import_r1_binding,
    import_r1_spec,
    numeric_close,
)
from fepdmc_toolkit.diagram_compiler.ir import KIND_ELECTRON, KIND_PHONON
from fepdmc_toolkit.scheme_d import exec as dexec
from fepdmc_toolkit.scheme_d.physics_first.contract import HolsteinRing, TwoBandHolstein
from fepdmc_toolkit.scheme_d.physics_first.index_oracle import (
    group_sum_matrix_twoband,
    group_sum_scalar,
)
from fepdmc_toolkit.scheme_d.physics_first.symbolic.bindings import (
    apply_scenario,
    base_binding,
)
from fepdmc_toolkit.scheme_d.physics_first.symbolic.plan_spec import scalar_spec, twoband_spec
from fepdmc_toolkit.scheme_d.physics_first.symbolic.r1 import r1_symbolic_build
from fepdmc_toolkit.scheme_d.physics_first.symbolic.r1_twoband import (
    r1_twoband_build,
    structural_eye,
)
from fepdmc_toolkit.scheme_d.physics_first.symbolic.reference_binder import ReferenceBinder

SCENARIOS = (
    "change_k",
    "change_q",
    "change_tau",
    "change_t",
    "change_omega",
    "change_g",
    "g_zero",
    "joint_parameters",
    "joint_sample",
    "repeat_q_then_separate",
    "repeat_dtau_then_separate",
    "change_L_analytic",
    "small_nonzero_g",
    "longer_legal_intervals",
    "same_input_repeat",
)


def _scalar_case(n: int, split: str, scenario: str):
    spec = scalar_spec(n)
    base = base_binding(spec, split)
    rt = base if scenario == "base" else apply_scenario(base, spec, scenario, 4242 + n)
    ir = import_r1_spec(spec)
    b = import_r1_binding(rt)
    plan = r1_symbolic_build(spec)
    env = ReferenceBinder(spec, plan.leaves).bind(rt)
    f_r1 = dexec.output_scalar(plan.graph, env, "F")
    ring = HolsteinRing(L=rt.L, t=rt.t, omega=rt.omega, g=rt.g)
    f_oracle = group_sum_scalar(n, rt.tau_array(), rt.q_array(), rt.k, ring)
    return ir, b, f_r1, f_oracle


@pytest.mark.parametrize("n", [1, 2, 3])
@pytest.mark.parametrize("split", ["DEV", "RESERVED_CONFIRMATION"])
@pytest.mark.parametrize("scenario", ["base", *SCENARIOS])
def test_scalar_evaluator_matches_r1_and_index_oracle(n, split, scenario):
    ir, b, f_r1, f_oracle = _scalar_case(n, split, scenario)
    f_mine = evaluate_diagramwise(ir, b)["F"]
    assert numeric_close(f_mine, f_r1)
    assert numeric_close(f_mine, f_oracle)


@pytest.mark.parametrize("n", [1, 2, 3])
def test_scalar_g_zero_is_exactly_zero(n):
    ir, b, f_r1, _ = _scalar_case(n, "DEV", "g_zero")
    f_mine = evaluate_diagramwise(ir, b)["F"]
    assert f_mine == 0.0 + 0j
    assert numeric_close(f_mine, f_r1)


def _twoband_case(n: int, split: str, scenario: str):
    spec = twoband_spec(n)
    base = base_binding(spec, split)
    rt = base if scenario == "base" else apply_scenario(base, spec, scenario, 9142 + n)
    ir = import_r1_spec(spec)
    b = import_r1_binding(rt)
    plan = r1_twoband_build(spec)
    env = ReferenceBinder(spec, plan.leaves).bind(rt)
    env.update(structural_eye())
    f_r1 = dexec.output_scalar(plan.graph, env, "F")
    op_r1 = dexec.interpret(plan.graph, env)["Op"]
    model = TwoBandHolstein(L=rt.L, t=rt.t, omega=rt.omega, g=rt.g, delta=rt.delta, gap=rt.gap)
    op_oracle = group_sum_matrix_twoband(n, rt.tau_array(), rt.q_array(), rt.k, model)
    return ir, b, f_r1, op_r1, op_oracle


@pytest.mark.parametrize("n", [1, 2])
@pytest.mark.parametrize("split", ["DEV", "RESERVED_CONFIRMATION"])
@pytest.mark.parametrize("scenario", ["base", *SCENARIOS])
def test_twoband_evaluator_matches_r1_and_index_oracle(n, split, scenario):
    ir, b, f_r1, op_r1, op_oracle = _twoband_case(n, split, scenario)
    out = evaluate_diagramwise(ir, b)
    assert numeric_close(out["F"], f_r1)
    assert numeric_close(out["Op"], op_r1)
    assert numeric_close(out["Op"], op_oracle)


def test_twoband_g_zero_op_is_exactly_zero():
    ir, b, f_r1, _, _ = _twoband_case(2, "DEV", "g_zero")
    out = evaluate_diagramwise(ir, b)
    assert np.all(out["Op"] == 0)
    assert numeric_close(out["F"], f_r1)


# ---------------------------------------------------------------------------
# Mutations: each fault must break agreement with the R1 oracle value.
# ---------------------------------------------------------------------------


def _oracle_value(n: int = 2, split: str = "DEV"):
    from dataclasses import replace

    spec = scalar_spec(n)
    rt = replace(
        base_binding(spec, split), L=8, k=np.pi / 4, q=(np.pi / 2, 3 * np.pi / 4), tau=(0.1, 0.27, 0.63, 0.94), g=0.8
    )
    plan = r1_symbolic_build(spec)
    f_r1 = dexec.output_scalar(plan.graph, ReferenceBinder(spec, plan.leaves).bind(rt), "F")
    return import_r1_spec(spec), import_r1_binding(rt), f_r1


def _rebind_diagram(ir: DiagramIR, index: int, diagram: Diagram) -> DiagramIR:
    diagrams = list(ir.diagrams)
    diagrams[index] = diagram
    return DiagramIR(ir.family, ir.variables, tuple(diagrams))


def test_mutation_sign_flip_breaks_agreement():
    ir, b, f_r1 = _oracle_value()
    d = ir.diagrams[0]
    mutated = Diagram(d.diagram_id, d.pairing, -1, d.multiplicity, d.factors)
    f_mut = evaluate_diagramwise(_rebind_diagram(ir, 0, mutated), b, check=False)["F"]
    assert not numeric_close(f_mut, f_r1)


def _flip_k_form(f: Factor) -> Factor:
    if f.kind != KIND_ELECTRON or f.k_form is None:
        return f
    if f.k_form.q_terms:
        k_form = MomentumForm(1, tuple((s, -c) for s, c in f.k_form.q_terms))
    else:
        k_form = MomentumForm(-1, ())
    return Factor(
        f.kind,
        interval=f.interval,
        open_chord_slots=f.open_chord_slots,
        k_form=k_form,
    )


def test_mutation_momentum_sign_breaks_agreement():
    ir, b, f_r1 = _oracle_value()
    d = ir.diagrams[0]
    mutated = Diagram(d.diagram_id, d.pairing, d.sign, d.multiplicity, tuple(_flip_k_form(f) for f in d.factors))
    f_mut = evaluate_diagramwise(_rebind_diagram(ir, 0, mutated), b, check=False)["F"]
    assert not numeric_close(f_mut, f_r1)


def test_mutation_phonon_interval_swap_breaks_agreement():
    ir, b, f_r1 = _oracle_value()
    d = ir.diagrams[0]
    factors = []
    flipped = False
    for f in d.factors:
        if f.kind == KIND_PHONON and not flipped:
            factors.append(
                Factor(f.kind, chord_slot=f.chord_slot, tau_open_idx=f.tau_close_idx, tau_close_idx=f.tau_open_idx)
            )
            flipped = True
        else:
            factors.append(f)
    mutated = Diagram(d.diagram_id, d.pairing, d.sign, d.multiplicity, tuple(factors))
    f_mut = evaluate_diagramwise(_rebind_diagram(ir, 0, mutated), b, check=False)["F"]
    assert not numeric_close(f_mut, f_r1)


def test_mutation_wrong_leaf_sharing_breaks_agreement():
    """Replacing one diagram's distinct electron factor by another diagram's
    leaf identity must change the value (wrong sharing is visible)."""
    ir, b, f_r1 = _oracle_value()
    d1, d2 = ir.diagrams[0], ir.diagrams[1]
    donor = next(f for f in d2.factors if f.kind == KIND_ELECTRON and f.interval == (1, 2))
    factors = tuple(donor if (f.kind == KIND_ELECTRON and f.interval == (1, 2)) else f for f in d1.factors)
    mutated = Diagram(d1.diagram_id, d1.pairing, d1.sign, d1.multiplicity, factors)
    f_mut = evaluate_diagramwise(_rebind_diagram(ir, 0, mutated), b, check=False)["F"]
    assert not numeric_close(f_mut, f_r1)


def test_mutation_prefactor_rejects_unknown_normalization():
    ir, b, _ = _oracle_value()
    from fepdmc_toolkit.diagram_compiler.evaluator import prefactor_value

    with pytest.raises(ValueError):
        prefactor_value(Factor("prefactor", normalization="g", exponent="2n"), b, ir.family.order)


def test_q_slot_permutation_changes_value_and_matches_r1():
    """Symbolic q-slot identity: swapping two distinct q values is a real
    physics change, not a numeric merge."""
    import dataclasses

    ir, b, f_r1 = _oracle_value()
    swapped = dataclasses.replace(b, q=(b.q[1], b.q[0]))
    f_mine = evaluate_diagramwise(ir, b)["F"]
    f_swap = evaluate_diagramwise(ir, swapped)["F"]
    assert numeric_close(f_mine, f_r1)
    assert not numeric_close(f_swap, f_r1)
    assert not numeric_close(f_swap, f_mine)


def test_no_scheme_d_import_outside_adapter():
    """The generated evaluator must not ride on R1 runtime helpers."""
    import inspect
    import sys

    import fepdmc_toolkit.diagram_compiler as pkg

    offenders = []
    for name, mod in list(sys.modules.items()):
        if name.startswith("fepdmc_toolkit.diagram_compiler") and mod is not None:
            src = inspect.getsource(sys.modules[name])
            for line in src.splitlines():
                s = line.strip()
                if s.startswith(("from fepdmc_toolkit.scheme_d", "import fepdmc_toolkit.scheme_d")):
                    offenders.append((name, s))
    assert not offenders, offenders
    assert pkg.__doc__ is not None


@pytest.mark.parametrize(
    "model,n",
    [("scalar", 1), ("scalar", 2), ("scalar", 3), ("scalar", 4), ("twoband", 1), ("twoband", 2), ("twoband", 3)],
)
def test_seeded_differential_all_lowerings(model, n):
    from fepdmc_toolkit.diagram_compiler import evaluate_dag, lower_diagramwise, lower_grouped
    from fepdmc_toolkit.scheme_d.physics_first.symbolic.guards import RuntimeBinding

    rng = np.random.default_rng(20260918 + n)
    spec = scalar_spec(n) if model == "scalar" else twoband_spec(n)
    ir = import_r1_spec(spec)
    dags = [lower_diagramwise(ir), lower_diagramwise(ir, share=True), lower_grouped(ir)]
    plan = r1_symbolic_build(spec) if model == "scalar" else r1_twoband_build(spec)
    binder = ReferenceBinder(spec, plan.leaves)
    for _ in range(12):
        L = int(rng.choice([3, 4, 8]))
        rt = RuntimeBinding(
            k=float(rng.integers(L) * 2 * np.pi / L),
            q=tuple(rng.integers(L, size=n) * 2 * np.pi / L),
            tau=tuple(np.cumsum(rng.uniform(0.01, 0.08, size=2 * n))),
            t=float(rng.uniform(0.3, 1.4)),
            omega=float(rng.uniform(0.4, 1.2)),
            g=float(rng.uniform(0.8, 1.6)),
            L=L,
            delta=float(rng.uniform(-0.6, 0.6)),
            gap=float(rng.uniform(0.2, 0.9)),
        )
        env = binder.bind(rt)
        if model == "twoband":
            env.update(structural_eye())
        reference = dexec.interpret(plan.graph, env)
        b = import_r1_binding(rt)
        naive = evaluate_diagramwise(ir, b)
        outputs = [naive] + [evaluate_dag(dag, b, ir) for dag in dags]
        if model == "scalar":
            physical = HolsteinRing(L=L, t=rt.t, omega=rt.omega, g=rt.g)
            independent = group_sum_scalar(n, rt.tau_array(), rt.q_array(), rt.k, physical)
            for out in outputs:
                assert numeric_close(out["F"], independent)
        else:
            physical = TwoBandHolstein(L=L, t=rt.t, omega=rt.omega, g=rt.g, delta=rt.delta, gap=rt.gap)
            independent = group_sum_matrix_twoband(n, rt.tau_array(), rt.q_array(), rt.k, physical)
            for out in outputs:
                assert numeric_close(out["Op"], independent)
        for out in outputs:
            assert numeric_close(out["F"], reference["F"])
            if model == "twoband":
                assert numeric_close(out["Op"], reference["Op"])


def test_twoband_wrong_multiplication_order_is_detected():
    from fepdmc_toolkit.scheme_d.physics_first.index_oracle import pairing_matrix_twoband, q_map_open_order

    ir, b, _, op_r1, _ = _twoband_case(2, "DEV", "base")
    model = TwoBandHolstein(L=b.L, t=b.t, omega=b.omega, g=b.g, delta=b.delta, gap=b.gap)
    wrong = sum(
        (
            pairing_matrix_twoband(
                d.pairing, b.tau_array(), q_map_open_order(d.pairing, b.q_array()), b.k, model, right_multiply_G=True
            )
            for d in ir.diagrams
        ),
        np.zeros((2, 2), dtype=complex),
    )
    assert not numeric_close(wrong, op_r1)
    assert numeric_close(evaluate_diagramwise(ir, b)["Op"], op_r1)
