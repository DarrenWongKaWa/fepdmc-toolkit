"""R1: open-phonon state DP. Returns DAG nodes. Does not import pairings."""

from __future__ import annotations

import time

import numpy as np

from ..ir import Binding, Graph, GraphBuilder, Layout, ObjectType
from .contract import HolsteinRing
from .meters import CostMeter


def r1_group_graph(
    n: int,
    taus: np.ndarray,
    qs: np.ndarray,
    k: float,
    ring: HolsteinRing,
) -> tuple[Graph, dict[str, np.ndarray], CostMeter]:
    if n < 1:
        raise ValueError("n>=1")
    if len(taus) != 2 * n or len(qs) != n:
        raise ValueError("x dimension")
    b = GraphBuilder(ObjectType.SHARED_X_GROUP.value, hash_cons=True)
    env: dict[str, np.ndarray] = {}
    meter = CostMeter(constructor="R1")
    memo: dict[tuple, int] = {}
    gv = ring.volume_vertex()
    pref = b.const(gv ** (2 * n), (1, 1), Layout.SCALAR.value)

    def g_leaf(k_eff: float, dtau: float) -> int:
        name = f"G_{k_eff:.16e}_{dtau:.16e}"
        if name not in env:
            env[name] = np.array([[np.exp(-ring.xi(k_eff) * dtau)]], dtype=np.complex128)
        return b.input(name, (1, 1), Layout.SCALAR.value, Binding.BOUND.value)

    def d_leaf(q: float, dtau: float) -> int:
        name = f"D_{q:.16e}_{dtau:.16e}"
        if name not in env:
            env[name] = np.array([[np.exp(-ring.omega * dtau)]], dtype=np.complex128)
        return b.input(name, (1, 1), Layout.SCALAR.value, Binding.BOUND.value)

    def rec(v: int, open_ph: tuple[tuple[int, int], ...], n_opened: int) -> int:
        key = (v, open_ph, n_opened)
        meter.recursive_calls += 1
        if key in memo:
            meter.memo_hits += 1
            return memo[key]
        if v == 2 * n + 1:
            node = b.const(1.0 if not open_ph else 0.0, (1, 1), Layout.SCALAR.value)
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
                dtau_ph = float(taus[v - 1] - taus[start - 1])
                dnode = d_leaf(float(qs[pid]), dtau_ph)
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
        k_eff = k - sum(float(qs[p]) for p, _ in open_after)
        dtau = float(taus[v] - taus[v - 1])
        return b.mul(g_leaf(k_eff, dtau), nxt)

    t0 = time.perf_counter()
    root = rec(1, (), 0)
    total = b.mul(pref, root)
    graph = b.finish({"F": total})
    meter.build_s = time.perf_counter() - t0
    meter.unique_states = len(memo)
    meter.dag_nodes = graph.n_nodes
    meter.dag_edges = sum(len(n.inputs) for n in graph.nodes)
    meter.leaf_bindings = len(env)
    meter.pairings_materialized = 0
    meter.notes = ("no_pairing_list",)
    return graph, env, meter
