"""Typed IR. Unique public writer: D03.

BOUND_C, SHARED_X_GROUP, and SUMMED_KERNEL are different mathematical
objects. A grouped or summed graph must not be fed to a bound-C weight
interface.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum


class ObjectType(str, Enum):
    BOUND_C = "BOUND_C"
    SHARED_X_GROUP = "SHARED_X_GROUP"
    SUMMED_KERNEL = "SUMMED_KERNEL"


class Layout(str, Enum):
    SCALAR = "scalar"
    DIAGONAL = "diagonal"
    DENSE = "dense"


class Binding(str, Enum):
    NONE = "none"
    FREE = "free"
    BOUND = "bound"
    SUMMED = "summed"


class Op(str, Enum):
    INPUT = "input"
    CONST = "const"
    ADD = "add"
    MUL = "mul"
    MATMUL = "matmul"
    SCALE_ROWS = "scale_rows"
    SCALE_COLS = "scale_cols"
    SUM = "sum"
    TRACE = "trace"
    TRANSPOSE = "transpose"
    CONJ = "conj"
    NEG = "neg"
    NORMALIZE = "normalize"
    LOG_NORM = "log_norm"


_OBJECT_TYPES = {e.value for e in ObjectType}


@dataclass(frozen=True)
class Node:
    op: str
    inputs: tuple[int, ...]
    layout: str
    shape: tuple[int, ...]
    axes: tuple[str, ...]
    binding: str
    var_name: str
    scale: complex
    phase: complex
    order: int
    resummed: bool
    material_id: str
    gauge_id: str
    attrs: tuple[tuple[str, str], ...]

    def plan_key(self) -> tuple:
        return (
            self.op,
            self.inputs,
            self.layout,
            self.shape,
            self.axes,
            self.binding,
            self.var_name,
            complex(self.scale),
            complex(self.phase),
            self.order,
            self.resummed,
            self.material_id,
            self.gauge_id,
            self.attrs,
        )


@dataclass(frozen=True)
class OpenPhonon:
    phonon_id: int
    q_name: str
    nu_name: str
    tau_start: str | None
    tau_end: str | None


@dataclass
class Graph:
    object_type: str
    nodes: list[Node]
    outputs: dict[str, int]
    bound_vars: tuple[str, ...]
    free_vars: tuple[str, ...]
    summed_vars: tuple[str, ...]
    open_phonons: tuple[OpenPhonon, ...]
    electron_in: str
    electron_out: str
    plan_id: str = ""
    notes: tuple[str, ...] = ()
    unit_inputs: tuple[str, ...] = ()

    @property
    def n_nodes(self) -> int:
        return len(self.nodes)

    def require_object(self, expected: str) -> None:
        if self.object_type != expected:
            raise TypeError(
                f"object_type={self.object_type} is not {expected}; "
                "BOUND_C, SHARED_X_GROUP, and SUMMED_KERNEL are not interchangeable"
            )


def as_bound_c_weight(graph: Graph) -> int:
    """Native D(C) hook. Grouped/summed graphs must raise."""
    graph.require_object(ObjectType.BOUND_C.value)
    if "weight" not in graph.outputs:
        raise KeyError("BOUND_C graph missing weight output")
    return graph.outputs["weight"]


class GraphBuilder:
    def __init__(
        self,
        object_type: str,
        *,
        hash_cons: bool = True,
        material_id: str = "synthetic",
        gauge_id: str = "default",
    ) -> None:
        if object_type not in _OBJECT_TYPES:
            raise ValueError(f"unknown object_type {object_type}")
        self.object_type = object_type
        self.hash_cons = hash_cons
        self.material_id = material_id
        self.gauge_id = gauge_id
        self.nodes: list[Node] = []
        self._cons: dict[tuple, int] = {}
        self.bound_vars: list[str] = []
        self.free_vars: list[str] = []
        self.summed_vars: list[str] = []
        self.open_phonons: list[OpenPhonon] = []
        self.electron_in = "k_in"
        self.electron_out = "k_out"
        self.outputs: dict[str, int] = {}
        self.notes: list[str] = []

    def add(self, node: Node) -> int:
        key = node.plan_key()
        if self.hash_cons and key in self._cons:
            return self._cons[key]
        idx = len(self.nodes)
        self.nodes.append(node)
        if self.hash_cons:
            self._cons[key] = idx
        return idx

    def _node(
        self,
        op: str,
        inputs: tuple[int, ...],
        layout: str,
        shape: tuple[int, ...],
        *,
        axes: tuple[str, ...] = (),
        binding: str = Binding.NONE.value,
        var_name: str = "",
        scale: complex = 1 + 0j,
        phase: complex = 1 + 0j,
        order: int = 0,
        resummed: bool = False,
        attrs: tuple[tuple[str, str], ...] = (),
    ) -> int:
        return self.add(
            Node(
                op=op,
                inputs=inputs,
                layout=layout,
                shape=shape,
                axes=axes,
                binding=binding,
                var_name=var_name,
                scale=complex(scale),
                phase=complex(phase),
                order=order,
                resummed=resummed,
                material_id=self.material_id,
                gauge_id=self.gauge_id,
                attrs=attrs,
            )
        )

    def input(
        self,
        var_name: str,
        shape: tuple[int, ...],
        layout: str,
        binding: str,
        *,
        axes: tuple[str, ...] = (),
        order: int = 0,
    ) -> int:
        if binding == Binding.BOUND.value and var_name not in self.bound_vars:
            self.bound_vars.append(var_name)
        elif binding == Binding.FREE.value and var_name not in self.free_vars:
            self.free_vars.append(var_name)
        elif binding == Binding.SUMMED.value and var_name not in self.summed_vars:
            self.summed_vars.append(var_name)
        return self._node(
            Op.INPUT.value,
            (),
            layout,
            shape,
            axes=axes,
            binding=binding,
            var_name=var_name,
            order=order,
        )

    def const(
        self,
        value: complex,
        shape: tuple[int, ...] = (1, 1),
        layout: str = Layout.SCALAR.value,
    ) -> int:
        return self._node(
            Op.CONST.value,
            (),
            layout,
            shape,
            attrs=(("value", repr(complex(value))),),
        )

    def add_op(self, a: int, b: int) -> int:
        na, nb = self.nodes[a], self.nodes[b]
        if na.shape != nb.shape:
            raise ValueError("add shape mismatch")
        layout = na.layout if na.layout == nb.layout else Layout.DENSE.value
        return self._node(Op.ADD.value, (a, b), layout, na.shape, axes=na.axes)

    def mul(self, a: int, b: int) -> int:
        na, nb = self.nodes[a], self.nodes[b]
        if na.shape != nb.shape:
            raise ValueError("mul shape mismatch")
        layout = na.layout if na.layout == nb.layout else Layout.DENSE.value
        return self._node(Op.MUL.value, (a, b), layout, na.shape, axes=na.axes)

    def matmul(self, a: int, b: int) -> int:
        na, nb = self.nodes[a], self.nodes[b]
        if len(na.shape) != 2 or len(nb.shape) != 2:
            raise ValueError("matmul expects rank-2")
        if na.shape[1] != nb.shape[0]:
            raise ValueError("matmul inner dimension mismatch")
        out_shape = (na.shape[0], nb.shape[1])
        if na.layout == Layout.SCALAR.value and nb.layout == Layout.SCALAR.value and out_shape == (1, 1):
            layout = Layout.SCALAR.value
        elif na.layout == Layout.DIAGONAL.value and nb.layout == Layout.DIAGONAL.value:
            layout = Layout.DIAGONAL.value
        else:
            layout = Layout.DENSE.value
        return self._node(Op.MATMUL.value, (a, b), layout, out_shape)

    def scale_rows(self, diag: int, mat: int) -> int:
        nd, nm = self.nodes[diag], self.nodes[mat]
        if nd.shape[0] != nm.shape[0]:
            raise ValueError("scale_rows mismatch")
        return self._node(Op.SCALE_ROWS.value, (diag, mat), nm.layout, nm.shape)

    def scale_cols(self, mat: int, diag: int) -> int:
        nm, nd = self.nodes[mat], self.nodes[diag]
        if nd.shape[0] != nm.shape[1]:
            raise ValueError("scale_cols mismatch")
        return self._node(Op.SCALE_COLS.value, (mat, diag), nm.layout, nm.shape)

    def trace(self, a: int) -> int:
        na = self.nodes[a]
        if len(na.shape) != 2 or na.shape[0] != na.shape[1]:
            raise ValueError("trace requires square")
        return self._node(Op.TRACE.value, (a,), Layout.SCALAR.value, (1, 1))

    def neg(self, a: int) -> int:
        na = self.nodes[a]
        return self._node(Op.NEG.value, (a,), na.layout, na.shape, axes=na.axes)

    def normalize(self, a: int) -> int:
        na = self.nodes[a]
        return self._node(Op.NORMALIZE.value, (a,), na.layout, na.shape, axes=na.axes)

    def log_norm(self, a: int) -> int:
        return self._node(Op.LOG_NORM.value, (a,), Layout.SCALAR.value, (1, 1))

    def finish(self, outputs: dict[str, int] | None = None) -> Graph:
        if outputs is not None:
            self.outputs = dict(outputs)
        return Graph(
            object_type=self.object_type,
            nodes=list(self.nodes),
            outputs=dict(self.outputs),
            bound_vars=tuple(self.bound_vars),
            free_vars=tuple(self.free_vars),
            summed_vars=tuple(self.summed_vars),
            open_phonons=tuple(self.open_phonons),
            electron_in=self.electron_in,
            electron_out=self.electron_out,
            plan_id=f"{self.object_type}:n{len(self.nodes)}",
            notes=tuple(self.notes),
        )


def layout_of_matrix(shape: tuple[int, ...], coupling: str) -> str:
    if shape == (1, 1) or coupling in {"constant_scalar", "momentum_scalar"}:
        if shape[0] == 1:
            return Layout.SCALAR.value
        if coupling != "noncommuting_matrix":
            return Layout.DIAGONAL.value
    if coupling == "noncommuting_matrix":
        return Layout.DENSE.value
    return Layout.DIAGONAL.value
