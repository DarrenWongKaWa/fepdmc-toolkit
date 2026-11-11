"""S08: symbolic post-enum / CSE / typed-CSE. Same leaf constructors as R1."""

from __future__ import annotations

import time

from fepdmc_toolkit.scheme_d.ir import GraphBuilder, Layout, ObjectType
from fepdmc_toolkit.scheme_d.physics_first.index_oracle import all_pairings
from fepdmc_toolkit.scheme_d.rewrites import RewriteRejected, fusion_consecutive_scales

from .cost import BuildReport, CallLedger
from .guards import refuse_unit_g_raw
from .leaves import LeafTable, electron_prop, graph_input, phonon_prop, prefactor_prop
from .plan_spec import PlanSpec
from .r1 import BuiltPlan


def _open_slots_after(pairing, v: int, pid_of: dict) -> tuple[int, ...]:
    slots = []
    for a, b in pairing:
        if a <= v < b:
            slots.append(pid_of[(a, b)])
    return tuple(sorted(slots))


def _pid_map(pairing) -> dict:
    ordered = tuple(sorted(pairing, key=lambda e: (e[0], e[1])))
    return {e: i for i, e in enumerate(ordered)}


def symbolic_enum_build(spec: PlanSpec, *, hash_cons: bool, ledger: CallLedger | None = None) -> BuiltPlan:
    if spec.model != "Holstein_scalar_vacuum":
        raise ValueError("scalar enum baseline")
    if ledger is not None:
        ledger.build_calls[spec.spec_id + ":" + ("CSE" if hash_cons else "enum")] += 1
    n = spec.n
    pairings = all_pairings(n)
    b = GraphBuilder(ObjectType.SHARED_X_GROUP.value, hash_cons=hash_cons)
    leaves = LeafTable()
    provider_d = spec.phonon_provider.replace("-", "")
    t0 = time.perf_counter()
    pref = graph_input(b, leaves, prefactor_prop(n=n, normalization_spec=spec.normalization_spec))
    terms = []
    for pairing in pairings:
        pid_of = _pid_map(pairing)
        node = b.const(1.0, (1, 1), Layout.SCALAR.value)
        verts = 2 * n
        for v in range(1, verts):
            gnode = graph_input(
                b,
                leaves,
                electron_prop(
                    provider="scalar_raw",
                    open_q_slots=_open_slots_after(pairing, v, pid_of),
                    tau_left_idx=v - 1,
                    tau_right_idx=v,
                    electron_interface=spec.electron_interface,
                ),
            )
            node = b.mul(gnode, node)
        for a, c in pairing:
            pid = pid_of[(a, c)]
            dnode = graph_input(
                b,
                leaves,
                phonon_prop(
                    provider=provider_d,
                    q_slot=pid,
                    mode_slot=0,
                    tau_open_idx=a - 1,
                    tau_close_idx=c - 1,
                ),
            )
            node = b.mul(dnode, node)
        terms.append(node)
    total = terms[0]
    for t in terms[1:]:
        total = b.add_op(total, t)
    graph = b.finish({"F": b.mul(pref, total)})
    ctor = "B-symbolic-CSE" if hash_cons else "B-symbolic-enum"
    report = BuildReport(
        spec_id=spec.spec_id,
        constructor=ctor,
        pairings_materialized=len(pairings),
        recursive_calls=0,
        unique_states=0,
        memo_hits=0,
        dag_nodes=graph.n_nodes,
        dag_edges=sum(len(n.inputs) for n in graph.nodes),
        n_leaves=len(leaves),
        open_boundary_width=0,
        build_s=time.perf_counter() - t0,
        notes=("enumerated_pairings", "symbolic_leaves"),
    )
    return BuiltPlan(spec=spec, graph=graph, leaves=leaves, report=report)


def typed_cse_or_reject_unit_g(plan: BuiltPlan) -> tuple[BuiltPlan, str]:
    """Ordinary CSE graph plus legal fusions. Unit-G skip is a failed guard, not an optimization."""
    fused = fusion_consecutive_scales(plan.graph)
    try:
        refuse_unit_g_raw(plan.spec)
        raise AssertionError("unit-G skip must not apply to raw physical_unshifted")
    except RewriteRejected:
        status = "G06_unit_G_rejected"
    report = BuildReport(
        spec_id=plan.spec.spec_id,
        constructor="B-symbolic-typed-CSE",
        pairings_materialized=plan.report.pairings_materialized,
        recursive_calls=0,
        unique_states=0,
        memo_hits=0,
        dag_nodes=fused.n_nodes,
        dag_edges=sum(len(n.inputs) for n in fused.nodes),
        n_leaves=len(plan.leaves),
        open_boundary_width=0,
        build_s=plan.report.build_s,
        notes=("typed_fusion_only", "unit_G_rejected", status),
    )
    return BuiltPlan(spec=plan.spec, graph=fused, leaves=plan.leaves, report=report), status
