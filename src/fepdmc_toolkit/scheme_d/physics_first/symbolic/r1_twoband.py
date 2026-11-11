"""Two-band R1: suffix operator R(v)=R(v+1) @ U @ V, matching the left-multiply oracle.

Untraced Op is an output. Trace-only agreement is not sufficient.
"""

from __future__ import annotations

import time

import numpy as np

from fepdmc_toolkit.scheme_d.ir import Binding, GraphBuilder, Layout, ObjectType

from .cost import BuildReport, CallLedger
from .leaves import LeafTable, electron_prop, graph_input, phonon_prop, vertex_prop
from .plan_spec import PlanSpec
from .r1 import BuiltPlan

STRUCT_EYE = "STRUCT_EYE_B2"


def structural_eye() -> dict[str, np.ndarray]:
    return {STRUCT_EYE: np.eye(2, dtype=np.complex128)}


def r1_twoband_build(spec: PlanSpec, ledger: CallLedger | None = None) -> BuiltPlan:
    if spec.model != "Holstein_twoband_hermitian":
        raise ValueError("two-band spec required")
    if spec.normalization_spec != "vertex_contains_g_over_sqrt_L":
        raise ValueError("do not add a second g/sqrt(L) prefactor")
    if ledger is not None:
        ledger.build_calls[spec.spec_id] += 1
    n = spec.n
    b = GraphBuilder(ObjectType.SHARED_X_GROUP.value, hash_cons=True)
    leaves = LeafTable()
    memo: dict[tuple, int] = {}
    recursive_calls = 0
    memo_hits = 0
    max_open = 0
    eye = b.input(STRUCT_EYE, (2, 2), Layout.DENSE.value, Binding.BOUND.value)
    zero = b.const(0.0, (2, 2), Layout.DENSE.value)
    provider_d = spec.phonon_provider.replace("-", "")

    def rec(v: int, open_ph: tuple[tuple[int, int], ...], n_opened: int) -> int:
        nonlocal recursive_calls, memo_hits, max_open
        key = (v, open_ph, n_opened)
        recursive_calls += 1
        if key in memo:
            memo_hits += 1
            return memo[key]
        max_open = max(max_open, len(open_ph))
        if v == 2 * n + 1:
            node = eye if not open_ph else zero
            memo[key] = node
            return node
        remaining = 2 * n - v + 1
        n_open = len(open_ph)
        terms: list[int] = []
        incoming = tuple(p for p, _ in open_ph)
        if n_opened < n and remaining - 1 >= n_open + 1:
            mnode = graph_input(
                b,
                leaves,
                vertex_prop(
                    provider=spec.vertex_provider,
                    direction="emit",
                    incoming_open_slots=incoming,
                    q_slot=n_opened,
                    mode_slot=0,
                    band_interface=spec.electron_interface,
                ),
            )
            terms.append(b.matmul(_after(v, open_ph + ((n_opened, v),), n_opened + 1), mnode))
        if n_open and remaining - 1 >= n_open - 1:
            for j, (pid, start) in enumerate(open_ph):
                new_open = open_ph[:j] + open_ph[j + 1 :]
                mnode = graph_input(
                    b,
                    leaves,
                    vertex_prop(
                        provider=spec.vertex_provider,
                        direction="absorb",
                        incoming_open_slots=incoming,
                        q_slot=pid,
                        mode_slot=0,
                        band_interface=spec.electron_interface,
                    ),
                )
                dnode = graph_input(
                    b,
                    leaves,
                    phonon_prop(
                        provider=provider_d,
                        q_slot=pid,
                        mode_slot=0,
                        tau_open_idx=start - 1,
                        tau_close_idx=v - 1,
                        layout_bands=2,
                    ),
                )
                inner = b.matmul(_after(v, new_open, n_opened), mnode)
                terms.append(b.matmul(dnode, inner))
        if not terms:
            node = zero
        else:
            node = terms[0]
            for t in terms[1:]:
                node = b.add_op(node, t)
        memo[key] = node
        return node

    def _after(v: int, open_after: tuple[tuple[int, int], ...], n_opened: int) -> int:
        nxt = rec(v + 1, open_after, n_opened)
        if v == 2 * n:
            return nxt
        unode = graph_input(
            b,
            leaves,
            electron_prop(
                provider="matrix_raw",
                open_q_slots=tuple(p for p, _ in open_after),
                tau_left_idx=v - 1,
                tau_right_idx=v,
                electron_interface=spec.electron_interface,
            ),
        )
        return b.matmul(nxt, unode)

    t0 = time.perf_counter()
    root = rec(1, (), 0)
    graph = b.finish({"Op": root, "F": b.trace(root)})
    report = BuildReport(
        spec_id=spec.spec_id,
        constructor="R1-symbolic-twoband",
        pairings_materialized=0,
        recursive_calls=recursive_calls,
        unique_states=len(memo),
        memo_hits=memo_hits,
        dag_nodes=graph.n_nodes,
        dag_edges=sum(len(n.inputs) for n in graph.nodes),
        n_leaves=len(leaves),
        open_boundary_width=max_open,
        build_s=time.perf_counter() - t0,
        notes=("no_pairing_list", "suffix_U_V", "no_extra_prefactor", "structural_eye", "untraced_Op"),
    )
    return BuiltPlan(spec=spec, graph=graph, leaves=leaves, report=report)
