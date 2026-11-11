"""Backend-ready NativeKernelIR: typed SSA values + basic-block CFG.

No Python closures, Binding objects, or exception types in the IR.
Stage 4 does not emit Fortran.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field

SCHEMA_VERSION = 1
IR_KIND = "NativeKernelIR"

TYPE_BOOL = "bool"
TYPE_I64 = "i64"
TYPE_F64 = "f64"
TYPE_C128 = "c128"
TYPE_C128_M2 = "c128_m2"
NATIVE_TYPES = (TYPE_BOOL, TYPE_I64, TYPE_F64, TYPE_C128, TYPE_C128_M2)

STATUS_OK = 0
INVALID_EXACT_TARGET = 1
INVALID_CHEAP_WEIGHT = 2
INVALID_PROPOSAL_RATIO = 3
INVALID_RANDOM_UNIFORM = 4
CACHE_KEY_MISMATCH = 5
INVALID_IR = 6
ERROR_NAMES = {
    STATUS_OK: "OK",
    INVALID_EXACT_TARGET: "INVALID_EXACT_TARGET",
    INVALID_CHEAP_WEIGHT: "INVALID_CHEAP_WEIGHT",
    INVALID_PROPOSAL_RATIO: "INVALID_PROPOSAL_RATIO",
    INVALID_RANDOM_UNIFORM: "INVALID_RANDOM_UNIFORM",
    CACHE_KEY_MISMATCH: "CACHE_KEY_MISMATCH",
    INVALID_IR: "INVALID_IR",
}

REJECT_NONE = 0
REJECT_STAGE1 = 1
REJECT_STAGE2 = 2

REGION_CHEAP = "cheap"
REGION_EXACT = "exact"
REGION_CONTROL = "control"


class NativeIRError(ValueError):
    def __init__(self, code: int, message: str) -> None:
        self.code = code
        super().__init__(f"{ERROR_NAMES.get(code, str(code))}: {message}")


@dataclass(frozen=True)
class LayoutField:
    name: str
    type: str
    shape: tuple[int, ...]

    def to_dict(self) -> dict:
        return {"name": self.name, "shape": list(self.shape), "type": self.type}

    @classmethod
    def from_dict(cls, payload: dict) -> LayoutField:
        return cls(payload["name"], payload["type"], tuple(payload["shape"]))


@dataclass(frozen=True)
class BindingLayout:
    schema_version: int
    model: str
    order: int
    fields: tuple[LayoutField, ...]

    def to_dict(self) -> dict:
        return {
            "fields": [f.to_dict() for f in self.fields],
            "model": self.model,
            "order": self.order,
            "schema_version": self.schema_version,
        }

    def to_json(self) -> str:
        return json.dumps(self.to_dict(), indent=2, sort_keys=True) + "\n"

    @classmethod
    def from_dict(cls, payload: dict) -> BindingLayout:
        if int(payload["schema_version"]) != SCHEMA_VERSION:
            raise NativeIRError(INVALID_IR, "unsupported binding-layout schema")
        fields = tuple(LayoutField.from_dict(item) for item in payload["fields"])
        return cls(SCHEMA_VERSION, str(payload["model"]), int(payload["order"]), fields)

    @classmethod
    def from_json(cls, text: str) -> BindingLayout:
        return cls.from_dict(json.loads(text))


def binding_layout_for(model: str, order: int) -> BindingLayout:
    fields = [
        LayoutField("k", TYPE_F64, ()),
        LayoutField("q", TYPE_F64, (order,)),
        LayoutField("tau", TYPE_F64, (2 * order,)),
        LayoutField("t", TYPE_F64, ()),
        LayoutField("omega", TYPE_F64, ()),
        LayoutField("g", TYPE_F64, ()),
        LayoutField("L", TYPE_I64, ()),
    ]
    if model == "Holstein_twoband_hermitian":
        fields.extend((LayoutField("delta", TYPE_F64, ()), LayoutField("gap", TYPE_F64, ())))
    return BindingLayout(SCHEMA_VERSION, model, order, tuple(fields))


@dataclass(frozen=True)
class NativeOp:
    id: str
    op: str
    type: str
    args: tuple[str, ...] = ()
    attrs: dict = field(default_factory=dict)
    region: str = REGION_CONTROL

    def to_dict(self) -> dict:
        payload = {"args": list(self.args), "id": self.id, "op": self.op, "region": self.region, "type": self.type}
        if self.attrs:
            payload["attrs"] = dict(self.attrs)
        return payload

    @classmethod
    def from_dict(cls, payload: dict) -> NativeOp:
        return cls(
            id=str(payload["id"]),
            op=str(payload["op"]),
            type=str(payload["type"]),
            args=tuple(payload.get("args", ())),
            attrs=dict(payload.get("attrs", {})),
            region=str(payload.get("region", REGION_CONTROL)),
        )


@dataclass(frozen=True)
class NativeTerminator:
    op: str
    args: tuple[str, ...] = ()
    attrs: dict = field(default_factory=dict)

    def to_dict(self) -> dict:
        payload = {"args": list(self.args), "op": self.op}
        if self.attrs:
            payload["attrs"] = dict(self.attrs)
        return payload

    @classmethod
    def from_dict(cls, payload: dict) -> NativeTerminator:
        return cls(str(payload["op"]), tuple(payload.get("args", ())), dict(payload.get("attrs", {})))


@dataclass(frozen=True)
class NativeBlock:
    id: str
    ops: tuple[NativeOp, ...]
    term: NativeTerminator

    def to_dict(self) -> dict:
        return {"id": self.id, "ops": [op.to_dict() for op in self.ops], "term": self.term.to_dict()}

    @classmethod
    def from_dict(cls, payload: dict) -> NativeBlock:
        ops = tuple(NativeOp.from_dict(item) for item in payload["ops"])
        return cls(str(payload["id"]), ops, NativeTerminator.from_dict(payload["term"]))


@dataclass(frozen=True)
class NativeGraph:
    id: str
    region: str
    ops: tuple[NativeOp, ...]
    output: str
    type: str

    def to_dict(self) -> dict:
        return {
            "id": self.id,
            "ops": [op.to_dict() for op in self.ops],
            "output": self.output,
            "region": self.region,
            "type": self.type,
        }

    @classmethod
    def from_dict(cls, payload: dict) -> NativeGraph:
        ops = tuple(NativeOp.from_dict(item) for item in payload["ops"])
        return cls(str(payload["id"]), str(payload["region"]), ops, str(payload["output"]), str(payload["type"]))


@dataclass(frozen=True)
class NativeKernelIR:
    schema_version: int
    ir_kind: str
    kind: str
    source: dict
    layout: BindingLayout
    graphs: tuple[NativeGraph, ...]
    blocks: tuple[NativeBlock, ...]
    entry: str

    def to_dict(self) -> dict:
        return {
            "blocks": [block.to_dict() for block in self.blocks],
            "entry": self.entry,
            "graphs": [graph.to_dict() for graph in self.graphs],
            "ir_kind": self.ir_kind,
            "kind": self.kind,
            "layout": self.layout.to_dict(),
            "schema_version": self.schema_version,
            "source": dict(self.source),
        }

    def to_json(self) -> str:
        return json.dumps(self.to_dict(), indent=2, sort_keys=True) + "\n"

    @classmethod
    def from_dict(cls, payload: dict) -> NativeKernelIR:
        if not isinstance(payload, dict) or payload.get("ir_kind") != IR_KIND:
            raise NativeIRError(INVALID_IR, "not a NativeKernelIR payload")
        extra = set(payload) - {
            "blocks",
            "entry",
            "graphs",
            "ir_kind",
            "kind",
            "layout",
            "schema_version",
            "source",
        }
        if extra:
            raise NativeIRError(INVALID_IR, f"unknown fields {sorted(extra)}")
        if int(payload["schema_version"]) != SCHEMA_VERSION:
            raise NativeIRError(INVALID_IR, "unsupported NativeKernelIR schema")
        graphs = tuple(NativeGraph.from_dict(item) for item in payload["graphs"])
        blocks = tuple(NativeBlock.from_dict(item) for item in payload["blocks"])
        ir = cls(
            SCHEMA_VERSION,
            IR_KIND,
            str(payload["kind"]),
            dict(payload["source"]),
            BindingLayout.from_dict(payload["layout"]),
            graphs,
            blocks,
            str(payload["entry"]),
        )
        if json.dumps(payload, sort_keys=True, allow_nan=False) != json.dumps(ir.to_dict(), sort_keys=True):
            raise NativeIRError(INVALID_IR, "unknown fields or noncanonical field types")
        return ir

    @classmethod
    def from_json(cls, text: str) -> NativeKernelIR:
        return cls.from_dict(json.loads(text))

    def graph(self, name: str) -> NativeGraph:
        for graph in self.graphs:
            if graph.id == name:
                return graph
        raise NativeIRError(INVALID_IR, f"missing graph {name}")

    def block(self, name: str) -> NativeBlock:
        for block in self.blocks:
            if block.id == name:
                return block
        raise NativeIRError(INVALID_IR, f"missing block {name}")
