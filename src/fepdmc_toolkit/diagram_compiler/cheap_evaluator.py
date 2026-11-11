"""Cheap evaluator generated from the same DiagramIR as the exact path.

Sibling of compile_evaluator: snapshot IR, lower with an inspectable
CheapPolicy, execute a scalar DAG. Never calls the exact interpreter,
historical R1, two-band eigh, vertex matrices, matmul, or trace.
"""

from __future__ import annotations

import json
import math
from dataclasses import dataclass

import numpy as np

from .canonical import validate_ir
from .cheap_cost import GraphCost, graph_cost
from .cheap_policy import CheapPolicy, CheapPolicyError, resolve_policy
from .evaluator import Binding, validate_binding
from .ir import (
    KIND_ELECTRON,
    KIND_PHONON,
    KIND_PREFACTOR,
    KIND_VERTEX,
    DiagramIR,
    MomentumForm,
)

P1_CORRESPONDENCE = {
    "cheap_semantics": "STATE_WEIGHT",
    "classes": (
        "representation_mismatch",
        "proposal_metadata_missing",
        "compiler_policy_different",
    ),
    "p1_semantics": "TRANSITION_SCORE",
    "verdict": "NOT_IDENTICAL",
}


class CheapEvalError(ValueError):
    """Cheap evaluation refused a binding or an illegal DAG op."""


@dataclass(frozen=True)
class CheapNode:
    op: str
    inputs: tuple[int, ...] = ()
    leaf: tuple | None = None
    const: object = None


@dataclass(frozen=True)
class CheapDag:
    family_model: str
    order: int
    policy_name: str
    nodes: tuple[CheapNode, ...]
    outputs: dict[str, int]


class _Builder:
    def __init__(self) -> None:
        self.nodes: list[CheapNode] = []
        self.by_key: dict[tuple, int] = {}

    def emit(self, op: str, key: tuple, *, inputs: tuple[int, ...] = (), leaf=None, const=None) -> int:
        got = self.by_key.get(key)
        if got is not None:
            return got
        nid = len(self.nodes)
        self.nodes.append(CheapNode(op, inputs, leaf, const))
        self.by_key[key] = nid
        return nid

    def const(self, key: tuple, value) -> int:
        return self.emit("const", key, const=value)

    def input_leaf(self, leaf: tuple) -> int:
        return self.emit("input", ("input", leaf), leaf=leaf)

    def mul(self, a: int, b: int) -> int:
        return self.emit("mul", ("mul", a, b), inputs=(a, b))

    def add(self, a: int, b: int) -> int:
        return self.emit("add", ("add", a, b), inputs=(a, b))

    def add_all(self, terms: list[int]) -> int:
        node = terms[0]
        for term in terms[1:]:
            node = self.add(node, term)
        return node


def _leaf_from_factor(factor) -> tuple:
    if factor.kind == KIND_ELECTRON:
        return (KIND_ELECTRON, factor.interval, factor.open_chord_slots, factor.k_form)
    if factor.kind == KIND_PHONON:
        return (KIND_PHONON, factor.chord_slot, factor.tau_open_idx, factor.tau_close_idx)
    if factor.kind == KIND_PREFACTOR:
        return (KIND_PREFACTOR, factor.normalization, factor.exponent)
    raise CheapPolicyError(f"uncovered category factor.{factor.kind}")


def _xi(k_eff: float, t: float) -> float:
    """Frozen xi_k=2t(1-cos k). Independent copy so cheap never calls electron_value."""
    return 2.0 * t * (1.0 - np.cos(k_eff))


def cheap_leaf_value(leaf: tuple, binding: Binding, order: int) -> complex:
    kind = leaf[0]
    if kind == KIND_ELECTRON:
        _, interval, _slots, k_form = leaf
        k_eff = float(k_form.k_coeff * binding.k)
        for slot, coeff in k_form.q_terms:
            k_eff += coeff * float(binding.q[slot])
        dtau = float(binding.tau[interval[1]] - binding.tau[interval[0]])
        return complex(np.exp(-_xi(k_eff, binding.t) * dtau))
    if kind == KIND_PHONON:
        _, _slot, tau_open, tau_close = leaf
        dtau = float(binding.tau[tau_close] - binding.tau[tau_open])
        return complex(np.exp(-binding.omega * dtau))
    if kind == KIND_PREFACTOR:
        _, normalization, exponent = leaf
        if normalization != "g_over_sqrt_L" or exponent != "2n":
            raise CheapEvalError(f"unsupported prefactor {normalization}/{exponent}")
        scale = binding.g / np.sqrt(binding.L)
        return complex(scale ** (2 * order))
    raise CheapEvalError(f"refused leaf kind {kind}")


_IMPLEMENTED_DROP = {
    "factor.vertex",
    "factor.prefactor",
    "factor.electron_propagator",
    "factor.phonon_propagator",
    "family.global_factors",
    "diagram.sign",
    "diagram.multiplicity",
    "eval.matmul",
    "eval.trace",
}


def lower_cheap(ir: DiagramIR, policy: CheapPolicy) -> CheapDag:
    validate_ir(ir)
    policy.check_ir(ir)
    if policy.action("eval.matmul") == "keep" or policy.action("eval.trace") == "keep":
        raise CheapPolicyError("cheap lowering refuses KEEP of matmul/trace")
    for category in policy.present_categories(ir):
        if policy.action(category) == "drop" and category not in _IMPLEMENTED_DROP:
            raise CheapPolicyError(f"DROP of {category} is not implemented; refuse rather than ignore")
    builder = _Builder()
    one = builder.const(("const_scalar_1",), complex(1.0))
    global_action = policy.action("family.global_factors")
    if global_action == "drop":
        pref = one
    elif ir.family.global_factors:
        (pref_factor,) = ir.family.global_factors
        pref = builder.input_leaf(_leaf_from_factor(pref_factor))
    else:
        # Vertices that carried g/sqrt(L) were dropped; restore the scalar scale.
        pref = builder.input_leaf((KIND_PREFACTOR, "g_over_sqrt_L", "2n"))
    terms: list[int] = []
    for diagram in ir.diagrams:
        weight = one
        for factor in diagram.factors:
            category = f"factor.{factor.kind}"
            try:
                action = policy.action(category)
            except CheapPolicyError as exc:
                raise CheapPolicyError(f"uncovered category {category}") from exc
            if action == "refuse":
                raise CheapPolicyError(f"refuse category {category}")
            if action == "drop":
                continue
            if factor.kind == KIND_VERTEX:
                raise CheapPolicyError("vertex cannot be kept in the cheap scalar DAG")
            if factor.kind not in (KIND_ELECTRON, KIND_PHONON, KIND_PREFACTOR):
                raise CheapPolicyError(f"uncovered category {category}")
            weight = builder.mul(weight, builder.input_leaf(_leaf_from_factor(factor)))
        if policy.action("diagram.sign") != "drop" and diagram.sign != 1:
            weight = builder.mul(weight, builder.const(("sign", diagram.sign), complex(diagram.sign)))
        if policy.action("diagram.multiplicity") != "drop" and diagram.multiplicity != 1:
            weight = builder.mul(
                weight,
                builder.const(("mult", diagram.multiplicity), complex(diagram.multiplicity)),
            )
        terms.append(weight)
    total = builder.add_all(terms)
    outputs = {"F_hat": builder.mul(pref, total)}
    return CheapDag(
        family_model=ir.family.model,
        order=ir.family.order,
        policy_name=policy.name,
        nodes=tuple(builder.nodes),
        outputs=outputs,
    )


def evaluate_cheap_dag(dag: CheapDag, binding: Binding) -> complex:
    cache: dict[int, complex] = {}
    for index, node in enumerate(dag.nodes):
        if node.op == "input":
            value = cheap_leaf_value(node.leaf, binding, dag.order)
        elif node.op == "const":
            value = node.const
        elif node.op == "add":
            value = cache[node.inputs[0]] + cache[node.inputs[1]]
        elif node.op == "mul":
            value = cache[node.inputs[0]] * cache[node.inputs[1]]
        elif node.op in ("matmul", "trace"):
            raise CheapEvalError(f"refused expensive op {node.op}")
        else:
            raise CheapEvalError(f"refused op {node.op}")
        cache[index] = value
    return complex(cache[dag.outputs["F_hat"]])


def odd_clip(value: float, limit: float) -> float:
    if limit <= 0:
        raise CheapEvalError("clip limit must be positive")
    if value > limit:
        return float(limit)
    if value < -limit:
        return float(-limit)
    return float(value)


def _log_positive_real(weight: complex) -> float:
    if not np.isfinite(weight.real) or not np.isfinite(weight.imag):
        raise CheapEvalError("nonfinite cheap weight")
    if abs(weight.imag) > 1e-12 or weight.real <= 0.0:
        raise CheapEvalError("non-positive cheap weight")
    return float(math.log(weight.real))


class CheapEvaluator:
    def __init__(self, ir: DiagramIR, policy: CheapPolicy, dag: CheapDag, cost: GraphCost) -> None:
        self.ir = ir
        self.policy = policy
        self.dag = dag
        self.cost = cost

    def evaluate(self, binding: Binding) -> dict:
        validate_binding(self.ir, binding)
        weight = evaluate_cheap_dag(self.dag, binding)
        log_w = _log_positive_real(weight)
        return {"F_hat": weight, "log_W_hat": log_w}

    def log_weight(self, binding: Binding) -> float:
        return float(self.evaluate(binding)["log_W_hat"])

    def transition_score(self, current: Binding, proposed: Binding) -> float:
        return self.log_weight(proposed) - self.log_weight(current)

    def to_dict(self) -> dict:
        return {
            "dag": _dag_to_dict(self.dag),
            "ir_kind": "CheapEvaluatorIR",
            "policy": self.policy.to_dict(),
            "schema_version": 1,
            "source_diagram_ir": self.ir.to_dict(),
        }

    def to_json(self) -> str:
        return json.dumps(self.to_dict(), indent=2, sort_keys=True) + "\n"

    @classmethod
    def from_dict(cls, payload: dict) -> CheapEvaluator:
        if not isinstance(payload, dict) or payload.get("ir_kind") != "CheapEvaluatorIR":
            raise CheapEvalError("not a CheapEvaluatorIR payload")
        extra = set(payload) - {"dag", "ir_kind", "policy", "schema_version", "source_diagram_ir"}
        if extra:
            raise CheapEvalError(f"unknown fields {sorted(extra)}")
        if payload.get("schema_version") != 1:
            raise CheapEvalError("unsupported CheapEvaluatorIR schema")
        ir = DiagramIR.from_dict(payload["source_diagram_ir"])
        policy = CheapPolicy.from_dict(payload["policy"])
        dag = _dag_from_dict(payload["dag"])
        if dag.order != ir.family.order or dag.family_model != ir.family.model:
            raise CheapEvalError("CheapEvaluatorIR DAG/IR family mismatch")
        cost = graph_cost(dag, expensive_electron=False)
        return cls(ir, policy, dag, cost)

    @classmethod
    def from_json(cls, text: str) -> CheapEvaluator:
        return cls.from_dict(json.loads(text))


def compile_cheap(ir: DiagramIR, policy: CheapPolicy | str | None = "propagator_only_v1") -> CheapEvaluator:
    resolved = resolve_policy(policy)
    snapshot = DiagramIR.from_json(ir.to_json())
    dag = lower_cheap(snapshot, resolved)
    cost = graph_cost(dag, expensive_electron=False)
    return CheapEvaluator(snapshot, resolved, dag, cost)


class ExactEvaluator:
    """Object wrapper around the unchanged stage-1 compile_evaluator callable."""

    def __init__(self, ir: DiagramIR) -> None:
        from .evaluator import compile_evaluator

        self.ir = DiagramIR.from_json(ir.to_json())
        self._evaluate = compile_evaluator(self.ir)

    def evaluate(self, binding: Binding) -> dict:
        return self._evaluate(binding)


def compile_exact(ir: DiagramIR) -> ExactEvaluator:
    return ExactEvaluator(ir)


def _leaf_to_dict(leaf: tuple) -> dict:
    kind = leaf[0]
    if kind == KIND_ELECTRON:
        _, interval, slots, k_form = leaf
        return {
            "interval": list(interval),
            "k_form": k_form.to_dict(),
            "kind": kind,
            "open_chord_slots": list(slots),
        }
    if kind == KIND_PHONON:
        _, slot, tau_open, tau_close = leaf
        return {
            "chord_slot": slot,
            "kind": kind,
            "tau_close_idx": tau_close,
            "tau_open_idx": tau_open,
        }
    if kind == KIND_PREFACTOR:
        _, normalization, exponent = leaf
        return {"exponent": exponent, "kind": kind, "normalization": normalization}
    raise CheapEvalError(f"refused leaf kind {kind}")


def _leaf_from_dict(payload: dict) -> tuple:
    kind = payload["kind"]
    if kind == KIND_ELECTRON:
        return (
            kind,
            tuple(payload["interval"]),
            tuple(payload["open_chord_slots"]),
            MomentumForm.from_dict(payload["k_form"]),
        )
    if kind == KIND_PHONON:
        return (kind, payload["chord_slot"], payload["tau_open_idx"], payload["tau_close_idx"])
    if kind == KIND_PREFACTOR:
        return (kind, payload["normalization"], payload["exponent"])
    raise CheapEvalError(f"refused leaf kind {kind}")


def _dag_to_dict(dag: CheapDag) -> dict:
    nodes = []
    for node in dag.nodes:
        item = {"inputs": list(node.inputs), "op": node.op}
        if node.leaf is not None:
            item["leaf"] = _leaf_to_dict(node.leaf)
        if node.const is not None:
            value = complex(node.const)
            item["const"] = [value.real, value.imag]
        nodes.append(item)
    return {
        "family_model": dag.family_model,
        "nodes": nodes,
        "order": dag.order,
        "outputs": dict(dag.outputs),
        "policy_name": dag.policy_name,
    }


def _dag_from_dict(payload: dict) -> CheapDag:
    nodes = []
    for item in payload["nodes"]:
        const = None
        if "const" in item:
            real, imag = item["const"]
            const = complex(real, imag)
        leaf = _leaf_from_dict(item["leaf"]) if "leaf" in item else None
        nodes.append(CheapNode(item["op"], tuple(item["inputs"]), leaf, const))
    return CheapDag(
        family_model=payload["family_model"],
        order=payload["order"],
        policy_name=payload["policy_name"],
        nodes=tuple(nodes),
        outputs=dict(payload["outputs"]),
    )
