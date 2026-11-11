"""DAG equivalence, recovered R1 node counts, sharing, and runtime isolation."""

from __future__ import annotations

import numpy as np
import pytest

from fepdmc_toolkit.diagram_compiler import (
    DagNode,
    EvalDag,
    evaluate_dag,
    evaluate_diagramwise,
    import_r1_binding,
    import_r1_spec,
    lower_diagramwise,
    lower_grouped,
    numeric_close,
)
from fepdmc_toolkit.diagram_compiler.ir import MomentumForm
from fepdmc_toolkit.scheme_d import exec as dexec
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


def _case(model: str, n: int, scenario: str = "base"):
    if model == "scalar":
        spec = scalar_spec(n)
    else:
        spec = twoband_spec(n)
    base = base_binding(spec, "DEV")
    rt = base if scenario == "base" else apply_scenario(base, spec, scenario, 777 + n)
    ir = import_r1_spec(spec)
    return ir, import_r1_binding(rt), spec, rt


@pytest.mark.parametrize("model,n", [("scalar", 1), ("scalar", 2), ("scalar", 3), ("twoband", 1), ("twoband", 2)])
def test_tree_dag_matches_naive_bitexact(model, n):
    ir, b, _, _ = _case(model, n)
    naive = evaluate_diagramwise(ir, b)
    dag = evaluate_dag(lower_diagramwise(ir, share=False), b, ir)
    assert dag["F"] == naive["F"]
    if model == "twoband":
        assert np.array_equal(dag["Op"], naive["Op"])


@pytest.mark.parametrize("model,n", [("scalar", 1), ("scalar", 2), ("scalar", 3), ("twoband", 1), ("twoband", 2)])
def test_shared_diagramwise_dag_matches_naive_bitexact(model, n):
    ir, b, _, _ = _case(model, n)
    naive = evaluate_diagramwise(ir, b)
    dag = evaluate_dag(lower_diagramwise(ir, share=True), b, ir)
    assert dag["F"] == naive["F"]
    if model == "twoband":
        assert np.array_equal(dag["Op"], naive["Op"])


@pytest.mark.parametrize(
    "model,n,scenario",
    [
        ("scalar", 2, "change_tau"),
        ("scalar", 3, "joint_parameters"),
        ("twoband", 2, "change_q"),
    ],
)
def test_grouped_dag_matches_naive_within_repo_tolerance(model, n, scenario):
    ir, b, _, _ = _case(model, n, scenario)
    naive = evaluate_diagramwise(ir, b)
    grouped = evaluate_dag(lower_grouped(ir), b, ir)
    assert numeric_close(grouped["F"], naive["F"])
    assert numeric_close(grouped["Op"], naive["Op"]) if model == "twoband" else True


@pytest.mark.parametrize("n", [1, 2, 3])
def test_grouped_scalar_node_and_leaf_counts_match_r1(n):
    ir, _, spec, _ = _case("scalar", n)
    plan = r1_symbolic_build(spec)
    stats = lower_grouped(ir).stats
    assert stats.unique_nodes == plan.report.dag_nodes
    assert stats.unique_leaves == plan.report.n_leaves


@pytest.mark.parametrize("n", [1, 2])
def test_grouped_twoband_node_and_leaf_counts_match_r1(n):
    ir, _, spec, _ = _case("twoband", n)
    plan = r1_twoband_build(spec)
    stats = lower_grouped(ir).stats
    assert stats.unique_nodes == plan.report.dag_nodes
    assert stats.unique_leaves == plan.report.n_leaves


@pytest.mark.parametrize("model,n", [("scalar", 2), ("scalar", 3), ("twoband", 2)])
def test_sharing_is_real(model, n):
    ir, _, _, _ = _case(model, n)
    tree = lower_diagramwise(ir, share=False)
    grouped = lower_grouped(ir)
    assert grouped.stats.unique_nodes < tree.stats.unique_nodes
    assert grouped.stats.leaf_uses > grouped.stats.unique_leaves
    assert grouped.stats.max_reuse >= 2


def test_lowering_is_deterministic():
    ir, _, _, _ = _case("scalar", 3)
    a = lower_grouped(ir)
    b = lower_grouped(ir)
    assert [nd.key for nd in a.nodes] == [nd.key for nd in b.nodes]
    assert a.outputs == b.outputs
    assert a.stats == b.stats


def test_no_state_leakage_between_bindings():
    ir, b1, _, rt = _case("scalar", 2)
    b2 = import_r1_binding(apply_scenario(rt, scalar_spec(2), "change_g", 31))
    dag = lower_grouped(ir)
    solo1 = evaluate_dag(dag, b1, ir)["F"]
    solo2 = evaluate_dag(dag, b2, ir)["F"]
    inter_a = evaluate_dag(dag, b1, ir)["F"]
    inter_b = evaluate_dag(dag, b2, ir)["F"]
    assert inter_a == solo1
    assert inter_b == solo2
    assert solo1 != solo2


def test_tampered_shared_leaf_identity_changes_value():
    """Swapping a shared leaf's identity in the grouped DAG must move the
    value: the shared structure is load-bearing, not decorative."""
    from dataclasses import replace

    ir, b, _, _ = _case("scalar", 2)
    b = replace(b, L=8, k=np.pi / 4, q=(np.pi / 2, 3 * np.pi / 4), tau=(0.1, 0.27, 0.63, 0.94), g=0.8)
    dag = lower_grouped(ir)
    correct = evaluate_dag(dag, b, ir)["F"]
    nodes = []
    patched = False
    for nd in dag.nodes:
        if (
            not patched
            and nd.op == "input"
            and nd.leaf[0] == "electron_propagator"
            and nd.leaf[1] == (1, 2)
            and len(nd.leaf[2]) == 2
        ):
            new_leaf = ("electron_propagator", (1, 2), (), MomentumForm(1, ()))
            nodes.append(DagNode(nd.op, ("input", new_leaf), nd.inputs, leaf=new_leaf, const=nd.const))
            patched = True
        else:
            nodes.append(nd)
    assert patched
    tampered = EvalDag(
        object_type=dag.object_type,
        family_model=dag.family_model,
        order=dag.order,
        share=dag.share,
        nodes=tuple(nodes),
        outputs=dict(dag.outputs),
        stats=dag.stats,
    )
    wrong = evaluate_dag(tampered, b, ir)["F"]
    assert not numeric_close(wrong, correct)


def test_compile_evaluator_snapshots_ir_and_lowers_once(monkeypatch):
    import fepdmc_toolkit.diagram_compiler.graph as graph
    from fepdmc_toolkit.diagram_compiler import compile_evaluator

    ir, a, spec, rt = _case("scalar", 2)
    b = import_r1_binding(apply_scenario(rt, spec, "change_g", 31))
    expected_a = evaluate_diagramwise(ir, a)
    expected_b = evaluate_diagramwise(ir, b)
    evaluator = compile_evaluator(ir)
    ir.diagrams[0].sign = -1

    def forbidden(*args, **kwargs):
        raise AssertionError("evaluation rebuilt a graph")

    monkeypatch.setattr(graph, "lower_diagramwise", forbidden)
    monkeypatch.setattr(graph, "lower_grouped", forbidden)
    assert evaluator(a) == expected_a
    assert evaluator(b) == expected_b
    assert evaluator(a) == expected_a


def test_dag_rejects_foreign_family():
    from fepdmc_toolkit.diagram_compiler import build_twoband_ir

    ir, b, _, _ = _case("scalar", 2)
    with pytest.raises(ValueError, match="DAG/IR"):
        evaluate_dag(lower_diagramwise(ir), b, build_twoband_ir(2))


def test_public_naive_evaluator_rejects_invalid_ir():
    ir, b, _, _ = _case("scalar", 2)
    ir.diagrams[0].sign = -1
    with pytest.raises(ValueError):
        evaluate_diagramwise(ir, b)


def test_vertex_leaf_key_requires_incoming_scope():
    from fepdmc_toolkit.diagram_compiler import Factor
    from fepdmc_toolkit.diagram_compiler.graph import leaf_key

    factor = Factor("vertex", vertex=1, chord_slot=0, direction="emit")
    with pytest.raises(ValueError):
        leaf_key(factor)
    assert leaf_key(factor, incoming=()) == ("vertex", "emit", (), 0)


@pytest.mark.parametrize("model", ["scalar", "twoband"])
def test_serialized_evaluator_runs_with_old_runtime_imports_blocked(tmp_path, model):
    import json
    import subprocess
    import sys
    from dataclasses import asdict
    from pathlib import Path

    ir, binding, spec, rt = _case(model, 2)
    plan = r1_symbolic_build(spec) if model == "scalar" else r1_twoband_build(spec)
    env = ReferenceBinder(spec, plan.leaves).bind(rt)
    if model == "twoband":
        env.update(structural_eye())
    reference = dexec.interpret(plan.graph, env)
    ir_path = tmp_path / "ir.json"
    binding_path = tmp_path / "binding.json"
    ir_path.write_text(ir.to_json())
    binding_path.write_text(json.dumps(asdict(binding)))
    program = """
import importlib.abc, json, sys
from pathlib import Path
class BlockOldRuntime(importlib.abc.MetaPathFinder):
    def find_spec(self, fullname, path=None, target=None):
        if fullname.startswith("fepdmc_toolkit.") and not fullname.startswith("fepdmc_toolkit.diagram_compiler"):
            raise AssertionError("old runtime import forbidden: " + fullname)
sys.meta_path.insert(0, BlockOldRuntime())
sys.path.insert(0, sys.argv[1])
from fepdmc_toolkit.diagram_compiler import Binding, DiagramIR, compile_evaluator
ir = DiagramIR.from_json(Path(sys.argv[2]).read_text())
out = compile_evaluator(ir)(Binding(**json.loads(Path(sys.argv[3]).read_text())))
print(json.dumps({key: {"real": value.real.tolist(), "imag": value.imag.tolist()}
                  if hasattr(value, "tolist") else {"real": value.real, "imag": value.imag}
                  for key, value in out.items()}))
"""
    src = Path(__file__).resolve().parents[1] / "src"
    result = subprocess.run(
        [sys.executable, "-I", "-c", program, str(src), str(ir_path), str(binding_path)],
        check=True,
        text=True,
        capture_output=True,
    )
    for key, value in json.loads(result.stdout).items():
        revived = np.asarray(value["real"]) + 1j * np.asarray(value["imag"])
        assert numeric_close(revived, reference[key])
