"""Evaluation semantics: binding, primitives, DAG interpreter.

The primitives here are written from the family contract in ir.py/canonical.py
(xi_k = 2t(1-cos k), Einstein vacuum phonon, g/sqrt(L) vertex). They do NOT
import or call the R1 evaluator (leaves.eval_leaf); equivalence to R1 is
established only by the differential tests.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .canonical import validate_ir
from .graph import EvalDag
from .ir import (
    KIND_ELECTRON,
    KIND_PHONON,
    KIND_PREFACTOR,
    KIND_VERTEX,
    MODEL_SCALAR,
    MODEL_TWOBAND,
    DiagramIR,
    Factor,
)

_ATOL = 1e-12
_RTOL = 1e-10
_TORUS_ATOL = 1e-12


class BindingError(ValueError):
    def __init__(self, code: str, message: str) -> None:
        self.code = code
        super().__init__(f"{code}: {message}")


@dataclass(frozen=True)
class Binding:
    """Runtime sample: one x point. Scope is declared by the IR variables."""

    k: float
    q: tuple[float, ...]
    tau: tuple[float, ...]
    t: float
    omega: float
    g: float
    L: int
    delta: float = 0.35
    gap: float = 0.5
    finite_model: bool = True
    continuous_extension: bool = False
    material_id: str = "synthetic"
    gauge_id: str = "default"
    data_version: str = "analytic_holstein_v1"

    def q_array(self) -> np.ndarray:
        return np.asarray(self.q, dtype=np.float64)

    def tau_array(self) -> np.ndarray:
        return np.asarray(self.tau, dtype=np.float64)


def _finite(*vals: float) -> bool:
    return all(np.isfinite(v) for v in vals)


def is_on_torus(value: float, L: int, atol: float = _TORUS_ATOL) -> bool:
    idx = int(np.round(value * L / (2.0 * np.pi)))
    snapped = 2.0 * np.pi * (idx % L) / L
    wrapped = 2.0 * np.pi * ((idx % L) + L) / L
    return bool(min(abs(value - snapped), abs(value - snapped - 2.0 * np.pi), abs(value - wrapped)) <= atol)


def torus_point(index: int, L: int) -> float:
    return float(2.0 * np.pi * (int(index) % L) / L)


def validate_binding(ir: DiagramIR, binding: Binding) -> None:
    n = ir.family.order
    q = binding.q_array()
    tau = binding.tau_array()
    if isinstance(binding.L, bool) or not isinstance(binding.L, (int, np.integer)) or int(binding.L) <= 0:
        raise BindingError("BC04", f"invalid_positive_L L={binding.L}")
    if (binding.material_id, binding.gauge_id, binding.data_version) != (
        "synthetic",
        "default",
        "analytic_holstein_v1",
    ):
        raise BindingError("BC08", "unsupported material, gauge, or data version")
    if (
        type(binding.finite_model) is not bool
        or type(binding.continuous_extension) is not bool
        or binding.finite_model == binding.continuous_extension
    ):
        raise BindingError("BC05", "choose exactly one finite-model or continuous-extension mode")
    if q.shape != (n,):
        raise BindingError("BC02", f"wrong_q_slot_count got {q.shape} want ({n},)")
    if tau.shape != (2 * n,):
        raise BindingError("BC02", f"wrong_tau_count got {tau.shape} want ({2 * n},)")
    fields = [binding.k, binding.t, binding.omega, binding.g, float(binding.L)]
    fields.extend(float(x) for x in q)
    fields.extend(float(x) for x in tau)
    if ir.family.model == MODEL_TWOBAND:
        fields.extend((binding.delta, binding.gap))
    if not _finite(*fields):
        raise BindingError("BC03", "nonfinite_input")
    if np.any(np.diff(tau) <= 0):
        raise BindingError("BC01", "unsorted_times; do not silently sort")
    if binding.finite_model:
        if not is_on_torus(binding.k, int(binding.L)):
            raise BindingError("BC05", "offgrid_k_as_finite_model")
        for i, qi in enumerate(q):
            if not is_on_torus(float(qi), int(binding.L)):
                raise BindingError("BC05", f"offgrid_q_as_finite_model slot={i}")


def numeric_close(value: complex | np.ndarray, reference: complex | np.ndarray) -> bool:
    """The repository tolerance standard: atol 1e-12 / rtol 1e-10."""
    v = np.asarray(value)
    r = np.asarray(reference)
    if not np.isfinite(v).all() or not np.isfinite(r).all():
        return False
    return bool(np.all(np.abs(v - r) <= _ATOL + _RTOL * np.abs(r)))


def _xi(k_eff: float, t: float) -> float:
    return 2.0 * t * (1.0 - np.cos(k_eff))


def electron_value(factor: Factor, binding: Binding, model: str) -> complex | np.ndarray:
    k_eff = float(factor.k_form.k_coeff * binding.k)
    for slot, coeff in factor.k_form.q_terms:
        k_eff += coeff * float(binding.q[slot])
    dtau = float(binding.tau[factor.interval[1]] - binding.tau[factor.interval[0]])
    if model == MODEL_SCALAR:
        return complex(np.exp(-_xi(k_eff, binding.t) * dtau))
    if model == MODEL_TWOBAND:
        xi = _xi(k_eff, binding.t)
        H = np.array(
            [[xi, binding.delta], [binding.delta, xi + binding.gap]],
            dtype=np.complex128,
        )
        w, v = np.linalg.eigh(H)
        return (v * np.exp(-w * dtau)) @ v.conj().T
    raise ValueError(model)


def phonon_value(factor: Factor, binding: Binding) -> complex:
    """Semantic value: scalar Einstein-vacuum D. The D*I2 leaf matrix used
    by R1's band-layout is an evaluation artifact, not the IR semantics."""
    dtau = float(binding.tau[factor.tau_close_idx] - binding.tau[factor.tau_open_idx])
    return complex(np.exp(-binding.omega * dtau))


def vertex_value(factor: Factor, binding: Binding, model: str) -> complex | np.ndarray:
    if model != MODEL_TWOBAND:
        raise ValueError("vertex factors exist only in the two-band family")
    gv = binding.g / np.sqrt(binding.L)
    return gv * np.array([[1.0, 0.0], [0.0, -1.0]], dtype=np.complex128)


def prefactor_value(factor: Factor, binding: Binding, order: int) -> complex:
    if factor.normalization != "g_over_sqrt_L":
        raise ValueError(factor.normalization)
    if factor.exponent != "2n":
        raise ValueError(factor.exponent)
    gv = binding.g / np.sqrt(binding.L)
    return complex(gv ** (2 * order))


def factor_value(factor: Factor, binding: Binding, ir: DiagramIR) -> complex | np.ndarray:
    kind = factor.kind
    if kind == KIND_ELECTRON:
        return electron_value(factor, binding, ir.family.model)
    if kind == KIND_PHONON:
        return phonon_value(factor, binding)
    if kind == KIND_VERTEX:
        return vertex_value(factor, binding, ir.family.model)
    if kind == KIND_PREFACTOR:
        return prefactor_value(factor, binding, ir.family.order)
    raise ValueError(kind)


def evaluate_diagramwise(ir: DiagramIR, binding: Binding, *, check: bool = True) -> dict:
    """Naive generated evaluator: no DAG, one product per diagram.

    Canonical semantics (documented in docs/diagram-compiler.md):
      scalar   F = prefactor * sum_C sign_C * mult_C * prod(factors of C)
      twoband  Op = sum_C (prod phonon scalars) * (vertex/electron factors
               left-multiplied in listed factor order), F = trace(Op)
    """
    if check:
        validate_ir(ir)
        validate_binding(ir, binding)
    model = ir.family.model
    n = ir.family.order
    if model == MODEL_SCALAR:
        (pref_f,) = ir.family.global_factors
        pref = prefactor_value(pref_f, binding, n)
        total = 0.0 + 0j
        for diagram in ir.diagrams:
            w = 1.0 + 0j
            for f in diagram.factors:
                w = w * factor_value(f, binding, ir)
            total += diagram.sign * diagram.multiplicity * w
        return {"F": complex(pref * total)}
    if model == MODEL_TWOBAND:
        op_total = np.zeros((2, 2), dtype=np.complex128)
        for diagram in ir.diagrams:
            state = np.eye(2, dtype=np.complex128)
            scale = 1.0 + 0j
            for f in diagram.factors:
                if f.kind == KIND_VERTEX:
                    state = vertex_value(f, binding, model) @ state
                elif f.kind == KIND_ELECTRON:
                    state = electron_value(f, binding, model) @ state
                elif f.kind == KIND_PHONON:
                    scale = scale * phonon_value(f, binding)
                else:
                    raise ValueError(f.kind)
            op_total += diagram.sign * diagram.multiplicity * (scale * state)
        return {"Op": op_total, "F": complex(np.trace(op_total))}
    raise ValueError(model)


def compile_evaluator(ir: DiagramIR, *, share: bool = True):
    """Lower an isolated IR snapshot once into an executable DAG closure."""
    from .graph import lower_diagramwise

    snapshot = DiagramIR.from_json(ir.to_json())
    dag = lower_diagramwise(snapshot, share=share)

    def evaluate(binding: Binding) -> dict:
        validate_binding(snapshot, binding)
        return evaluate_dag(dag, binding, snapshot, check=False)

    return evaluate


def evaluate_dag(dag: EvalDag, binding: Binding, ir: DiagramIR, *, check: bool = True) -> dict:
    """Execute an internally lowered DAG, computing shared nodes once."""
    if check:
        if (dag.object_type, dag.family_model, dag.order) != (ir.object_type, ir.family.model, ir.family.order):
            raise ValueError("DAG/IR family mismatch")
        validate_ir(ir)
        validate_binding(ir, binding)
    cache: dict[int, complex | np.ndarray] = {}
    for i, node in enumerate(dag.nodes):
        if node.op == "input":
            value = factor_value(_factor_from_leaf(node.leaf, ir), binding, ir)
        elif node.op == "const":
            value = node.const
        elif node.op == "add":
            value = cache[node.inputs[0]] + cache[node.inputs[1]]
        elif node.op == "mul":
            value = cache[node.inputs[0]] * cache[node.inputs[1]]
        elif node.op == "matmul":
            a = cache[node.inputs[0]]
            b = cache[node.inputs[1]]
            if np.ndim(a) == 0 or np.ndim(b) == 0:
                # Phonon D is a scalar multiple of the identity; the grouped
                # lowering keeps R1's matmul node shape with scalar values.
                value = a * b
            else:
                value = a @ b
        elif node.op == "trace":
            value = complex(np.trace(cache[node.inputs[0]]))
        else:
            raise ValueError(node.op)
        cache[i] = value
    return {name: cache[i] for name, i in dag.outputs.items()}


def _factor_from_leaf(leaf: tuple, ir: DiagramIR) -> Factor:
    kind = leaf[0]
    if kind == KIND_ELECTRON:
        _, interval, slots, k_form = leaf
        return Factor(kind, interval=interval, open_chord_slots=slots, k_form=k_form)
    if kind == KIND_PHONON:
        _, slot, tau_o, tau_c = leaf
        return Factor(kind, chord_slot=slot, tau_open_idx=tau_o, tau_close_idx=tau_c)
    if kind == KIND_VERTEX:
        # Leaf identity tuple: (kind, direction, incoming_slots, chord_slot).
        _, direction, _incoming, slot = leaf
        return Factor(kind, chord_slot=slot, direction=direction)
    if kind == KIND_PREFACTOR:
        _, normalization, exponent = leaf
        return Factor(kind, normalization=normalization, exponent=exponent)
    raise ValueError(kind)
