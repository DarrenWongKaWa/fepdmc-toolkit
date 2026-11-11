"""Explicit-factor DAG lowering and a separate recovered R1 recurrence."""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass

import numpy as np

from .canonical import (
    KIND_ELECTRON,
    KIND_PHONON,
    KIND_PREFACTOR,
    KIND_VERTEX,
    chord_slots,
    validate_ir,
)
from .ir import (
    MODEL_SCALAR,
    MODEL_TWOBAND,
    DiagramIR,
    Factor,
    MomentumForm,
)


@dataclass(frozen=True, eq=False)
class DagNode:
    op: str
    key: tuple
    inputs: tuple[int, ...] = ()
    leaf: tuple | None = None
    const: object = None


@dataclass
class DagStats:
    unique_nodes: int
    node_uses: int
    unique_leaves: int
    leaf_uses: int
    max_reuse: int


@dataclass
class EvalDag:
    object_type: str
    family_model: str
    order: int
    share: bool
    nodes: tuple[DagNode, ...]
    outputs: dict[str, int]
    stats: DagStats


def leaf_key(f: Factor, *, incoming: tuple[int, ...] | None = None) -> tuple:
    kind = f.kind
    if kind == KIND_ELECTRON:
        return (kind, f.interval, f.open_chord_slots, f.k_form)
    if kind == KIND_PHONON:
        return (kind, f.chord_slot, f.tau_open_idx, f.tau_close_idx)
    if kind == KIND_VERTEX:
        if incoming is None:
            raise ValueError("vertex leaf needs incoming open slots")
        return (kind, f.direction, incoming, f.chord_slot)
    if kind == KIND_PREFACTOR:
        return (kind, f.normalization, f.exponent)
    raise ValueError(kind)


class _Builder:
    """Ordered structural CSE; identity constants remain shared in either mode."""

    def __init__(self, share: bool) -> None:
        self.share = share
        self.nodes: list[DagNode] = []
        self.by_key: dict[tuple, int] = {}

    def emit(self, node: DagNode) -> int:
        if self.share:
            got = self.by_key.get(node.key)
            if got is not None:
                return got
        nid = len(self.nodes)
        self.nodes.append(node)
        if self.share:
            self.by_key[node.key] = nid
        return nid

    def input_leaf(self, f: Factor, *, incoming: tuple[int, ...] | None = None) -> int:
        key = leaf_key(f, incoming=incoming)
        return self.emit(DagNode("input", ("input", key), leaf=key))

    def const_node(self, key: tuple, value) -> int:
        return self.emit(DagNode("const", key, const=value))

    def add(self, a: int, b: int) -> int:
        return self.emit(DagNode("add", ("add", a, b), inputs=(a, b)))

    def mul(self, a: int, b: int) -> int:
        return self.emit(DagNode("mul", ("mul", a, b), inputs=(a, b)))

    def matmul(self, a: int, b: int) -> int:
        return self.emit(DagNode("matmul", ("matmul", a, b), inputs=(a, b)))

    def trace(self, a: int) -> int:
        return self.emit(DagNode("trace", ("trace", a), inputs=(a,)))

    def add_all(self, terms: list[int]) -> int:
        node = terms[0]
        for t in terms[1:]:
            node = self.add(node, t)
        return node

    def finish(self, outputs: dict[str, int], *, object_type: str, model: str, order: int, share: bool) -> EvalDag:
        fanout: defaultdict[int, int] = defaultdict(int)
        node_uses = 0
        for nd in self.nodes:
            for i in nd.inputs:
                fanout[i] += 1
                node_uses += 1
        unique_leaves = sum(1 for nd in self.nodes if nd.op == "input")
        leaf_uses = sum(f for i, f in fanout.items() if self.nodes[i].op == "input")
        stats = DagStats(
            unique_nodes=len(self.nodes),
            node_uses=node_uses,
            unique_leaves=unique_leaves,
            leaf_uses=leaf_uses,
            max_reuse=max(fanout.values(), default=0),
        )
        return EvalDag(
            object_type=object_type,
            family_model=model,
            order=order,
            share=share,
            nodes=tuple(self.nodes),
            outputs=dict(outputs),
            stats=stats,
        )


def _electron_factor(v: int, slots: tuple[int, ...]) -> Factor:
    return Factor(
        KIND_ELECTRON,
        interval=(v - 1, v),
        open_chord_slots=slots,
        k_form=MomentumForm(1, tuple((s, -1) for s in slots)),
    )


def _phonon_factor(slot: int, start: int, close: int) -> Factor:
    return Factor(KIND_PHONON, chord_slot=slot, tau_open_idx=start - 1, tau_close_idx=close - 1)


def _vertex_factor(v: int, slot: int, direction: str) -> Factor:
    return Factor(KIND_VERTEX, vertex=v, chord_slot=slot, direction=direction)


def lower_diagramwise(ir: DiagramIR, *, share: bool = False) -> EvalDag:
    """Naive semantics as a DAG: per-diagram chain, summed in IR order."""
    validate_ir(ir)
    b = _Builder(share)
    model = ir.family.model
    n = ir.family.order
    if model == MODEL_SCALAR:
        one = b.const_node(("const_scalar_1",), complex(1.0))
        terms = []
        for diagram in ir.diagrams:
            w = one
            for f in diagram.factors:
                w = b.mul(w, b.input_leaf(f))
            terms.append(w)
        total = b.add_all(terms)
        (pref_f,) = ir.family.global_factors
        pref = b.input_leaf(pref_f)
        outputs = {"F": b.mul(pref, total)}
    elif model == MODEL_TWOBAND:
        eye = b.const_node(("const_eye2",), np.eye(2, dtype=np.complex128))
        one = b.const_node(("const_scalar_1",), complex(1.0))
        terms = []
        for diagram in ir.diagrams:
            state = eye
            scale = one
            slots = chord_slots(diagram.pairing)
            incoming_by_v = {
                v: tuple(sorted(slots[e] for e in diagram.pairing if e[0] < v <= e[1])) for v in range(1, 2 * n + 1)
            }
            for f in diagram.factors:
                if f.kind == KIND_VERTEX:
                    node = b.input_leaf(f, incoming=incoming_by_v[f.vertex])
                    state = b.matmul(node, state)
                elif f.kind == KIND_ELECTRON:
                    node = b.input_leaf(f)
                    state = b.matmul(node, state)
                elif f.kind == KIND_PHONON:
                    scale = b.mul(scale, b.input_leaf(f))
                else:
                    raise ValueError(f.kind)
            terms.append(b.mul(scale, state))
        op = b.add_all(terms)
        outputs = {"Op": op, "F": b.trace(op)}
    else:
        raise ValueError(model)
    return b.finish(outputs, object_type=ir.object_type, model=model, order=n, share=share)


def lower_grouped(ir: DiagramIR) -> EvalDag:
    """Recovered R1 representation: DP over (v, open_ph, n_opened).

    Mirrors the implicit R1 recursion branch-for-branch (open/close
    sufficiency conditions, slot-ordered closes, left-assoc sums,
    gnode*suffix / suffix@U_v folding) so the DAG structure is directly
    comparable with the historical R1 Graph.
    """
    validate_ir(ir)
    b = _Builder(share=True)
    model = ir.family.model
    n = ir.family.order
    memo: dict[tuple, int] = {}

    def remove(chord: tuple[int, int], open_ph: tuple[tuple[int, int], ...]) -> tuple[tuple[int, int], ...]:
        j = open_ph.index(chord)
        return open_ph[:j] + open_ph[j + 1 :]

    if model == MODEL_SCALAR:
        pref_f = ir.family.global_factors[0]
        if pref_f.kind != KIND_PREFACTOR:
            raise ValueError("scalar family must carry a prefactor")
        pref = b.input_leaf(pref_f)
        # R1 emits scalar consts lazily inside the recursion; mirror that so
        # node-count comparisons are not skewed by unused pre-emitted consts.
        one_key, zero_key = ("const_scalar_1",), ("const_scalar_0",)
        one_val, zero_val = complex(1.0), complex(0.0)

        def rec(v: int, open_ph: tuple[tuple[int, int], ...], n_opened: int) -> int:
            key = (v, open_ph, n_opened)
            if key in memo:
                return memo[key]
            if v == 2 * n + 1:
                ok = (not open_ph) and n_opened == n
                node = b.const_node(one_key, one_val) if ok else b.const_node(zero_key, zero_val)
                memo[key] = node
                return node
            remaining = 2 * n - v + 1
            n_open = len(open_ph)
            terms: list[int] = []
            if n_opened < n and remaining - 1 >= n_open + 1:
                terms.append(_after(v, open_ph + ((n_opened, v),), n_opened + 1))
            if n_open and remaining - 1 >= n_open - 1:
                for chord in open_ph:
                    pid, start = chord
                    dnode = b.input_leaf(_phonon_factor(pid, start, v))
                    terms.append(b.mul(dnode, _after(v, remove(chord, open_ph), n_opened)))
            if not terms:
                node = b.const_node(zero_key, zero_val)
            else:
                node = b.add_all(terms)
            memo[key] = node
            return node

        def _after(v: int, open_after: tuple[tuple[int, int], ...], n_opened: int) -> int:
            nxt = rec(v + 1, open_after, n_opened)
            if v == 2 * n:
                return nxt
            slots = tuple(p for p, _ in open_after)
            gnode = b.input_leaf(_electron_factor(v, slots))
            return b.mul(gnode, nxt)

        root = rec(1, (), 0)
        outputs = {"F": b.mul(pref, root)}
        return b.finish(outputs, object_type=ir.object_type, model=model, order=n, share=True)

    if model == MODEL_TWOBAND:
        eye = b.const_node(("const_eye2",), np.eye(2, dtype=np.complex128))
        zero = b.const_node(("const_zero2",), np.zeros((2, 2), dtype=np.complex128))

        def rec(v: int, open_ph: tuple[tuple[int, int], ...], n_opened: int) -> int:
            key = (v, open_ph, n_opened)
            if key in memo:
                return memo[key]
            if v == 2 * n + 1:
                node = eye if not open_ph else zero
                memo[key] = node
                return node
            remaining = 2 * n - v + 1
            n_open = len(open_ph)
            incoming = tuple(p for p, _ in open_ph)
            terms: list[int] = []
            if n_opened < n and remaining - 1 >= n_open + 1:
                mnode = b.input_leaf(_vertex_factor(v, n_opened, "emit"), incoming=incoming)
                terms.append(b.matmul(_after(v, open_ph + ((n_opened, v),), n_opened + 1), mnode))
            if n_open and remaining - 1 >= n_open - 1:
                for chord in open_ph:
                    pid, start = chord
                    mnode = b.input_leaf(_vertex_factor(v, pid, "absorb"), incoming=incoming)
                    dnode = b.input_leaf(_phonon_factor(pid, start, v))
                    inner = b.matmul(_after(v, remove(chord, open_ph), n_opened), mnode)
                    terms.append(b.matmul(dnode, inner))
            node = zero if not terms else b.add_all(terms)
            memo[key] = node
            return node

        def _after(v: int, open_after: tuple[tuple[int, int], ...], n_opened: int) -> int:
            nxt = rec(v + 1, open_after, n_opened)
            if v == 2 * n:
                return nxt
            slots = tuple(p for p, _ in open_after)
            unode = b.input_leaf(_electron_factor(v, slots))
            return b.matmul(nxt, unode)

        root = rec(1, (), 0)
        outputs = {"Op": root, "F": b.trace(root)}
        return b.finish(outputs, object_type=ir.object_type, model=model, order=n, share=True)

    raise ValueError(model)
