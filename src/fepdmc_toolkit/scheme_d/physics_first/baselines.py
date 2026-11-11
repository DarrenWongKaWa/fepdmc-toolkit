"""Five baselines. Enumeration lives here, not inside R1–R4."""

from __future__ import annotations

import time

import numpy as np

from ..exec import output_scalar
from ..ir import Binding, GraphBuilder, Layout, ObjectType
from .contract import HolsteinRing
from .index_oracle import all_pairings, is_proper, pairing_weight_scalar, q_map_open_order
from .meters import CostMeter


def _g_name(k_eff, dtau):
    return f"G_{k_eff:.16e}_{dtau:.16e}"


def enum_sum_graph(
    n: int,
    taus: np.ndarray,
    qs: np.ndarray,
    k: float,
    ring: HolsteinRing,
    *,
    hash_cons: bool,
    proper_only: bool = False,
) -> tuple:
    b = GraphBuilder(ObjectType.SHARED_X_GROUP.value, hash_cons=hash_cons)
    env = {}
    meter = CostMeter(constructor="B-enum-CSE" if hash_cons else "B-enum")
    pairings = all_pairings(n)
    if proper_only:
        pairings = tuple(p for p in pairings if is_proper(p))
    meter.pairings_materialized = len(pairings)
    pref = b.const(ring.volume_vertex() ** (2 * n), (1, 1), Layout.SCALAR.value)

    def g_leaf(k_eff, dtau):
        name = _g_name(k_eff, dtau)
        if name not in env:
            env[name] = np.array([[np.exp(-ring.xi(k_eff) * dtau)]], dtype=np.complex128)
        return b.input(name, (1, 1), Layout.SCALAR.value, Binding.BOUND.value)

    def d_leaf(q, dtau):
        name = f"D_{q:.16e}_{dtau:.16e}"
        if name not in env:
            env[name] = np.array([[np.exp(-ring.omega * dtau)]], dtype=np.complex128)
        return b.input(name, (1, 1), Layout.SCALAR.value, Binding.BOUND.value)

    t0 = time.perf_counter()
    terms = []
    for pairing in pairings:
        q_of = q_map_open_order(pairing, qs)
        node = b.const(1.0, (1, 1), Layout.SCALAR.value)
        verts = 2 * n
        for v in range(1, verts):
            open_e = [e for e in pairing if e[0] <= v < e[1]]
            k_eff = k - sum(q_of[e] for e in open_e)
            dtau = float(taus[v] - taus[v - 1])
            node = b.mul(g_leaf(k_eff, dtau), node)
        for e in pairing:
            dtau_ph = float(taus[e[1] - 1] - taus[e[0] - 1])
            node = b.mul(d_leaf(q_of[e], dtau_ph), node)
        terms.append(node)
    total = terms[0]
    for t in terms[1:]:
        total = b.add_op(total, t)
    graph = b.finish({"F": b.mul(pref, total)})
    meter.build_s = time.perf_counter() - t0
    meter.dag_nodes = graph.n_nodes
    meter.leaf_bindings = len(env)
    return graph, env, meter


def b_enum_direct(n, taus, qs, k, ring) -> complex:
    acc = 0.0 + 0j
    for p in all_pairings(n):
        acc += pairing_weight_scalar(p, taus, q_map_open_order(p, qs), k, ring)
    return complex(acc)


def evaluate_graph(graph, env) -> complex:
    return output_scalar(graph, env, "F")
