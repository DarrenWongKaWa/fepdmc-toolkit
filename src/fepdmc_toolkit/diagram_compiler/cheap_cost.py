"""Static cheap-vs-exact cost model. Node counts, not wall time."""

from __future__ import annotations

from dataclasses import dataclass

from .ir import KIND_ELECTRON, KIND_PHONON, KIND_PREFACTOR, KIND_VERTEX, MODEL_TWOBAND


@dataclass(frozen=True)
class GraphCost:
    unique_nodes: int
    n_electron: int
    n_phonon: int
    n_vertex: int
    n_prefactor: int
    n_add: int
    n_mul: int
    n_matmul: int
    n_trace: int
    n_eigh: int
    expensive_primitives: int

    def to_dict(self) -> dict:
        return {
            "expensive_primitives": self.expensive_primitives,
            "n_add": self.n_add,
            "n_eigh": self.n_eigh,
            "n_electron": self.n_electron,
            "n_matmul": self.n_matmul,
            "n_mul": self.n_mul,
            "n_phonon": self.n_phonon,
            "n_prefactor": self.n_prefactor,
            "n_trace": self.n_trace,
            "n_vertex": self.n_vertex,
            "unique_nodes": self.unique_nodes,
        }


def _leaf_kind(node) -> str | None:
    leaf = getattr(node, "leaf", None)
    if leaf is None:
        return None
    if isinstance(leaf, tuple) and leaf:
        return leaf[0]
    if isinstance(leaf, dict):
        return leaf.get("kind")
    return None


def graph_cost(dag, *, expensive_electron: bool) -> GraphCost:
    n_electron = n_phonon = n_vertex = n_prefactor = 0
    n_add = n_mul = n_matmul = n_trace = n_eigh = 0
    nodes = dag.nodes
    for node in nodes:
        op = node.op
        kind = _leaf_kind(node)
        if op == "input":
            if kind == KIND_ELECTRON:
                n_electron += 1
                if expensive_electron:
                    n_eigh += 1
            elif kind == KIND_PHONON:
                n_phonon += 1
            elif kind == KIND_VERTEX:
                n_vertex += 1
            elif kind == KIND_PREFACTOR:
                n_prefactor += 1
        elif op == "add":
            n_add += 1
        elif op == "mul":
            n_mul += 1
        elif op == "matmul":
            n_matmul += 1
        elif op == "trace":
            n_trace += 1
        elif op not in ("const", "input"):
            raise ValueError(f"unrecognized DAG op {op}")
    expensive = n_eigh + n_vertex + n_matmul + n_trace
    return GraphCost(
        unique_nodes=len(nodes),
        n_electron=n_electron,
        n_phonon=n_phonon,
        n_vertex=n_vertex,
        n_prefactor=n_prefactor,
        n_add=n_add,
        n_mul=n_mul,
        n_matmul=n_matmul,
        n_trace=n_trace,
        n_eigh=n_eigh,
        expensive_primitives=expensive,
    )


def lower_exact_cost(ir) -> GraphCost:
    from .graph import lower_diagramwise

    dag = lower_diagramwise(ir, share=True)
    return graph_cost(dag, expensive_electron=ir.family.model == MODEL_TWOBAND)
