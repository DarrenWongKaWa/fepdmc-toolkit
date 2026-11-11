"""S06 unique ABI writer: preallocated packed buffer. Not dict-pack billed as zero-cost."""

from __future__ import annotations

import ctypes
import time

import numpy as np

from fepdmc_toolkit.scheme_d.fortran_ir import CompiledGraph

from .cost import CallLedger
from .guards import GuardError, RuntimeBinding, validate_binding
from .leaves import LeafTable, eval_leaf
from .plan_spec import PlanSpec


class EvalContext:
    """Own x/work/y buffers. Compiled library is read-only shared."""

    def __init__(self, compiled: CompiledGraph) -> None:
        self.compiled = compiled
        self.fn = compiled.fn
        self.x = np.empty(max(1, compiled.ir.n_input), dtype=np.complex128)
        self.work = np.empty(max(1, compiled.ir.n_work), dtype=np.complex128)
        self.output = np.empty(max(1, compiled.ir.n_output), dtype=np.complex128)
        self.status = ctypes.c_int()

    def run_filled(self) -> np.ndarray:
        self.fn(
            self.x.ctypes.data,
            self.work.ctypes.data,
            self.output.ctypes.data,
            ctypes.byref(self.status),
        )
        if self.status.value:
            raise FloatingPointError(f"generated evaluator status={self.status.value}")
        return self.output.copy()


class PackedBinder:
    def __init__(
        self,
        spec: PlanSpec,
        leaves: LeafTable,
        compiled: CompiledGraph,
        ledger: CallLedger | None = None,
        structural: dict[str, np.ndarray] | None = None,
    ) -> None:
        self.spec = spec
        self.leaves = leaves
        self.compiled = compiled
        self.ledger = ledger
        self.structural = structural or {}
        self.names = list(compiled.ir.variables)
        self.offsets = {n: compiled.ir.variables[n]["offset"] for n in self.names}
        self.sizes = {n: compiled.ir.variables[n]["size"] for n in self.names}
        missing = [n for n in self.names if n not in leaves.descriptors and n not in self.structural]
        if missing:
            raise KeyError(f"ABI names not in leaf table or structural: {missing}")

    def fill(self, binding: RuntimeBinding, ctx: EvalContext) -> float:
        validate_binding(self.spec, binding)
        t0 = time.perf_counter()
        if self.ledger is not None:
            self.ledger.packed_bind_calls[self.spec.spec_id] += 1
        x = ctx.x
        for name in self.names:
            item = self.compiled.ir.variables[name]
            off = item["offset"]
            size = item["size"]
            if name in self.structural:
                value = np.asarray(self.structural[name], dtype=np.complex128)
            else:
                value = eval_leaf(self.leaves.descriptors[name], binding)
            flat = np.asarray(value, dtype=np.complex128).ravel(order="F")
            if flat.size != size:
                raise ValueError(f"leaf {name} size {flat.size} != ABI {size}")
            if not np.isfinite(flat).all():
                raise GuardError("G03", f"nonfinite leaf {name}")
            x[off : off + size] = flat
        return time.perf_counter() - t0

    def evaluate(self, binding: RuntimeBinding, ctx: EvalContext) -> tuple[complex, float, float]:
        try:
            bind_s = self.fill(binding, ctx)
        except GuardError:
            if self.ledger is not None:
                self.ledger.reject_calls[self.spec.spec_id] += 1
            raise
        t1 = time.perf_counter()
        out = ctx.run_filled()
        eval_s = time.perf_counter() - t1
        if self.ledger is not None:
            self.ledger.eval_calls[self.spec.spec_id] += 1
        return complex(out[0]), bind_s, eval_s
