"""S03: R1 open-phonon DP on a PlanSpec. No sample numbers, no pairing list."""

from __future__ import annotations

import time
from dataclasses import dataclass

from fepdmc_toolkit.scheme_d.ir import Graph, GraphBuilder, Layout, ObjectType

from .cost import BuildReport, CallLedger
from .leaves import LeafTable, electron_prop, graph_input, phonon_prop, prefactor_prop
from .plan_spec import PlanSpec


@dataclass
class BuiltPlan:
    spec: PlanSpec
    graph: Graph
    leaves: LeafTable
    report: BuildReport

    def fingerprint(self) -> tuple:
        return tuple((n.op, n.inputs, n.var_name, n.layout, n.shape, n.attrs) for n in self.graph.nodes)


def r1_symbolic_build(spec: PlanSpec, ledger: CallLedger | None = None) -> BuiltPlan:
    if spec.model != "Holstein_scalar_vacuum":
        raise ValueError("r1_symbolic_build is scalar; use r1_twoband")
    if spec.propagator_convention != "physical_unshifted":
        raise ValueError(spec.propagator_convention)
    if spec.n < 1:
        raise ValueError("n>=1")
    if ledger is not None:
        ledger.build_calls[spec.spec_id] += 1
    n = spec.n
    b = GraphBuilder(ObjectType.SHARED_X_GROUP.value, hash_cons=True)
    leaves = LeafTable()
    memo: dict[tuple, int] = {}
    recursive_calls = 0
    memo_hits = 0
    max_open = 0
    provider_e = "scalar_raw"
    provider_d = spec.phonon_provider.replace("-", "")

    pref = graph_input(b, leaves, prefactor_prop(n=n, normalization_spec=spec.normalization_spec))

    # 1-based vertices v=1..2n. Handoff R(0) is rec(1); handoff v==2n is v==2n+1 here.
    def rec(v: int, open_ph: tuple[tuple[int, int], ...], n_opened: int) -> int:
        nonlocal recursive_calls, memo_hits, max_open
        key = (v, open_ph, n_opened)
        recursive_calls += 1
        if key in memo:
            memo_hits += 1
            return memo[key]
        max_open = max(max_open, len(open_ph))
        if v == 2 * n + 1:
            ok = (not open_ph) and n_opened == n
            node = b.const(1.0 if ok else 0.0, (1, 1), Layout.SCALAR.value)
            memo[key] = node
            return node
        remaining = 2 * n - v + 1
        n_open = len(open_ph)
        terms: list[int] = []
        if n_opened < n and remaining - 1 >= n_open + 1:
            terms.append(_after(v, open_ph + ((n_opened, v),), n_opened + 1))
        if n_open and remaining - 1 >= n_open - 1:
            for j, (pid, start) in enumerate(open_ph):
                new_open = open_ph[:j] + open_ph[j + 1 :]
                dnode = graph_input(
                    b,
                    leaves,
                    phonon_prop(
                        provider=provider_d,
                        q_slot=pid,
                        mode_slot=0,
                        tau_open_idx=start - 1,
                        tau_close_idx=v - 1,
                    ),
                )
                terms.append(b.mul(dnode, _after(v, new_open, n_opened)))
        if not terms:
            node = b.const(0.0, (1, 1), Layout.SCALAR.value)
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
        gnode = graph_input(
            b,
            leaves,
            electron_prop(
                provider=provider_e,
                open_q_slots=tuple(p for p, _ in open_after),
                tau_left_idx=v - 1,
                tau_right_idx=v,
                electron_interface=spec.electron_interface,
            ),
        )
        return b.mul(gnode, nxt)

    t0 = time.perf_counter()
    root = rec(1, (), 0)
    total = b.mul(pref, root)
    graph = b.finish({"F": total})
    build_s = time.perf_counter() - t0
    report = BuildReport(
        spec_id=spec.spec_id,
        constructor="R1-symbolic",
        pairings_materialized=0,
        recursive_calls=recursive_calls,
        unique_states=len(memo),
        memo_hits=memo_hits,
        dag_nodes=graph.n_nodes,
        dag_edges=sum(len(n.inputs) for n in graph.nodes),
        n_leaves=len(leaves),
        open_boundary_width=max_open,
        build_s=build_s,
        notes=("no_pairing_list", "symbolic_leaves", "spec_only_builder"),
    )
    return BuiltPlan(spec=spec, graph=graph, leaves=leaves, report=report)


class PlanCache:
    def __init__(self, ledger: CallLedger) -> None:
        self.ledger = ledger
        self._plans: dict[tuple, BuiltPlan] = {}

    def get(self, spec: PlanSpec, builder=r1_symbolic_build) -> BuiltPlan:
        key = spec.identity_key() + (builder.__name__,)
        cached = self._plans.get(key)
        if cached is not None:
            return cached
        plan = builder(spec, ledger=self.ledger)
        self._plans[key] = plan
        return plan
