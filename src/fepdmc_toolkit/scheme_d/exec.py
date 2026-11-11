"""IR interpreter. MATMUL is left-associative and non-commuting."""

from __future__ import annotations

from collections.abc import Mapping

import numpy as np

from .ir import Graph, Layout, Op


def _as_matrix(x: np.ndarray, shape: tuple[int, ...]) -> np.ndarray:
    a = np.asarray(x)
    if a.shape == shape:
        return a
    if a.shape == () and shape == (1, 1):
        return np.array([[complex(a)]], dtype=np.complex128)
    if a.ndim == 1 and shape == (a.shape[0], a.shape[0]):
        return np.diag(a).astype(np.complex128)
    raise ValueError(f"cannot reshape {a.shape} to {shape}")


def interpret(graph: Graph, env: Mapping[str, np.ndarray]) -> dict[str, np.ndarray]:
    validate_bindings(graph, env)
    vals: list[np.ndarray] = []
    for node in graph.nodes:
        if node.op == Op.INPUT.value:
            if node.var_name not in env:
                raise KeyError(f"missing input {node.var_name}")
            v = _as_matrix(env[node.var_name], node.shape)
        elif node.op == Op.CONST.value:
            raw = dict(node.attrs)["value"]
            v = np.full(node.shape, complex(raw), dtype=np.complex128)
        elif node.op == Op.ADD.value:
            v = vals[node.inputs[0]] + vals[node.inputs[1]]
        elif node.op == Op.MUL.value:
            v = vals[node.inputs[0]] * vals[node.inputs[1]]
        elif node.op == Op.MATMUL.value:
            left = vals[node.inputs[0]]
            right = vals[node.inputs[1]]
            v = left @ right
        elif node.op == Op.SCALE_ROWS.value:
            d = np.diag(vals[node.inputs[0]]) if vals[node.inputs[0]].ndim == 2 else vals[node.inputs[0]]
            v = np.asarray(d, dtype=np.complex128)[:, None] * vals[node.inputs[1]]
        elif node.op == Op.SCALE_COLS.value:
            d = np.diag(vals[node.inputs[1]]) if vals[node.inputs[1]].ndim == 2 else vals[node.inputs[1]]
            v = vals[node.inputs[0]] * np.asarray(d, dtype=np.complex128)[None, :]
        elif node.op == Op.TRACE.value:
            v = np.array([[np.trace(vals[node.inputs[0]])]], dtype=np.complex128)
        elif node.op == Op.NEG.value:
            v = -vals[node.inputs[0]]
        elif node.op == Op.TRANSPOSE.value:
            v = vals[node.inputs[0]].T
        elif node.op == Op.CONJ.value:
            v = np.conjugate(vals[node.inputs[0]])
        elif node.op in {Op.NORMALIZE.value, Op.LOG_NORM.value}:
            source = vals[node.inputs[0]]
            norm = float(np.max(np.abs(source)))
            if not np.isfinite(norm) or norm == 0:
                raise FloatingPointError("normalization requires a finite nonzero norm")
            v = source / norm if node.op == Op.NORMALIZE.value else np.array([[np.log(norm)]], dtype=np.complex128)
        else:
            raise ValueError(f"unknown op {node.op}")
        if node.scale != 1 + 0j:
            v = v * node.scale
        if node.phase != 1 + 0j:
            v = v * node.phase
        vals.append(np.asarray(v, dtype=np.complex128))
    return {name: vals[idx] for name, idx in graph.outputs.items()}


def validate_bindings(graph: Graph, env: Mapping[str, np.ndarray]) -> None:
    """Type facts used by exact rewrites are runtime contracts, not guesses."""
    for name in graph.unit_inputs:
        if name not in env or not np.array_equal(np.asarray(env[name]), np.ones((1, 1))):
            raise ValueError(f"unit-propagator guard failed: {name}")
    for node in graph.nodes:
        if node.op != Op.INPUT.value:
            continue
        a = _as_matrix(env[node.var_name], node.shape)
        if not np.isfinite(a).all():
            raise ValueError(f"nonfinite input: {node.var_name}")
        if node.layout == Layout.DIAGONAL.value and not np.array_equal(a, np.diag(np.diag(a))):
            raise ValueError(f"diagonal input contract failed: {node.var_name}")


def output_scalar(graph: Graph, env: Mapping[str, np.ndarray], name: str = "weight") -> complex:
    out = interpret(graph, env)[name]
    return complex(np.asarray(out).reshape(-1)[0])


def count_ops(graph: Graph) -> dict[str, int]:
    acc: dict[str, int] = {}
    for n in graph.nodes:
        acc[n.op] = acc.get(n.op, 0) + 1
    acc["n_nodes"] = graph.n_nodes
    return acc
