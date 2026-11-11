"""S07: one Fortran compile per PlanSpec identity. Multiple safe contexts."""

from __future__ import annotations

from pathlib import Path

from fepdmc_toolkit.scheme_d.fortran_ir import CompiledGraph
from fepdmc_toolkit.scheme_d.ir import Graph

from .cost import CallLedger
from .fused_kernel import FusedKernel
from .leaves import LeafTable
from .packed_binder import EvalContext
from .plan_spec import PlanSpec


class CompileCache:
    def __init__(self, root: Path, ledger: CallLedger) -> None:
        self.root = Path(root)
        self.root.mkdir(parents=True, exist_ok=True)
        self.ledger = ledger
        self._compiled: dict[tuple, CompiledGraph] = {}
        self._fused: dict[tuple, FusedKernel] = {}
        self._source: dict[tuple, str] = {}

    def _key(self, spec: PlanSpec, tag: str) -> tuple:
        return spec.identity_key() + (tag,)

    def compile(self, spec: PlanSpec, graph: Graph, *, tag: str = "R1-symbolic") -> CompiledGraph:
        key = self._key(spec, tag)
        hit = self._compiled.get(key)
        if hit is not None:
            return hit
        self.ledger.compile_calls[spec.spec_id + ":" + tag] += 1
        directory = self.root / f"{spec.spec_id}_{tag}"
        compiled = CompiledGraph(graph, directory)
        self._compiled[key] = compiled
        self._source[key] = compiled.ir.source
        return compiled

    def source(self, spec: PlanSpec, *, tag: str = "R1-symbolic") -> str:
        return self._source[self._key(spec, tag)]

    def spawn_context(self, spec: PlanSpec, *, tag: str = "R1-symbolic") -> EvalContext:
        return EvalContext(self._compiled[self._key(spec, tag)])

    def compile_fused(self, spec: PlanSpec, graph: Graph, leaves: LeafTable, *, tag: str) -> FusedKernel:
        key = self._key(spec, tag + ":fused")
        hit = self._fused.get(key)
        if hit is not None:
            return hit
        self.ledger.compile_calls[spec.spec_id + ":" + tag + ":fused"] += 1
        directory = self.root / f"{spec.spec_id}_{tag}_fused"
        kernel = FusedKernel(spec, graph, leaves, directory, tag=tag)
        self._fused[key] = kernel
        self._source[key] = kernel.source
        return kernel
