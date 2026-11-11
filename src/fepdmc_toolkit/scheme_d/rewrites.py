"""Exact rewrites with preserved scale/phase, metadata and binding guards."""

from __future__ import annotations

from dataclasses import replace

from .ir import Binding, Graph, GraphBuilder, Layout, Op


class RewriteRejected(ValueError):
    """Input is outside the exact applicability domain."""


def scalar_unit_propagator_applicable(shape, energy_shift):
    return shape in {(1,), (1, 1)} and complex(energy_shift) == 0j


def _copy_node(b, node, remap):
    return b.add(replace(node, inputs=tuple(remap[i] for i in node.inputs)))


def _finish(graph, b, remap, note, **kwargs):
    return replace(
        graph,
        nodes=list(b.nodes),
        outputs={k: remap[v] for k, v in graph.outputs.items()},
        notes=graph.notes + (note,),
        **kwargs,
    )


def apply_scalar_unit_skip(graph: Graph, energy_shifts: dict[str, complex]) -> Graph:
    guards = set(graph.unit_inputs)
    for node in graph.nodes:
        if node.op == Op.INPUT.value and node.var_name.startswith("G"):
            if not scalar_unit_propagator_applicable(node.shape, energy_shifts.get(node.var_name, 1)):
                raise RewriteRejected(f"scalar unit-G skip inapplicable for {node.var_name} shape={node.shape}")
            guards.add(node.var_name)
    b = GraphBuilder(graph.object_type, hash_cons=True)
    remap = {}
    for i, node in enumerate(graph.nodes):
        if node.op == Op.INPUT.value and node.var_name.startswith("G"):
            remap[i] = b.add(
                replace(
                    node,
                    op=Op.CONST.value,
                    inputs=(),
                    binding=Binding.NONE.value,
                    var_name="",
                    attrs=(("value", "(1+0j)"),),
                )
            )
        else:
            remap[i] = _copy_node(b, node, remap)
    return _finish(graph, b, remap, "scalar_unit_skip", unit_inputs=tuple(sorted(guards)))


def lower_diagonal_matmul(graph: Graph) -> Graph:
    b = GraphBuilder(graph.object_type, hash_cons=True)
    remap = {}
    for i, node in enumerate(graph.nodes):
        if node.op == Op.MATMUL.value:
            left, right = (graph.nodes[j] for j in node.inputs)
            if left.layout == Layout.DIAGONAL.value and right.layout != Layout.SCALAR.value:
                remap[i] = b.add(replace(node, op=Op.SCALE_ROWS.value, inputs=tuple(remap[j] for j in node.inputs)))
                continue
            if right.layout == Layout.DIAGONAL.value and left.layout != Layout.SCALAR.value:
                remap[i] = b.add(replace(node, op=Op.SCALE_COLS.value, inputs=tuple(remap[j] for j in node.inputs)))
                continue
        remap[i] = _copy_node(b, node, remap)
    return _finish(graph, b, remap, "diagonal_lowering")


def fusion_consecutive_scales(graph: Graph) -> Graph:
    b = GraphBuilder(graph.object_type, hash_cons=True)
    remap = {}
    fused = 0
    for i, node in enumerate(graph.nodes):
        if node.op == Op.SCALE_ROWS.value:
            inner = graph.nodes[node.inputs[1]]
            if inner.op == Op.SCALE_ROWS.value:
                dprod = b.mul(remap[node.inputs[0]], remap[inner.inputs[0]])
                remap[i] = b.add(
                    replace(
                        node,
                        inputs=(dprod, remap[inner.inputs[1]]),
                        scale=node.scale * inner.scale,
                        phase=node.phase * inner.phase,
                    )
                )
                fused += 1
                continue
        remap[i] = _copy_node(b, node, remap)
    return _finish(graph, b, remap, f"fused_scales={fused}")
