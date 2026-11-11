"""Public compile-once session. Create once; evaluate many times. Not per-sample build."""

from __future__ import annotations

import time
from pathlib import Path

import numpy as np

from .baselines import symbolic_enum_build
from .compile_cache import CompileCache
from .cost import CallLedger
from .fused_kernel import FusedContext, FusedKernel
from .guards import GuardError, RuntimeBinding, validate_binding
from .plan_spec import PlanSpec
from .r1 import BuiltPlan, PlanCache, r1_symbolic_build


def _cse_builder(spec: PlanSpec, ledger: CallLedger | None = None):
    return symbolic_enum_build(spec, hash_cons=True, ledger=ledger)


_BUILDERS = {
    "R1": r1_symbolic_build,
    "R1-symbolic": r1_symbolic_build,
    "CSE": _cse_builder,
    "B-symbolic-CSE": _cse_builder,
}

# Code cache only. Sessions are not interned: each has its own _ctx and guard contract.
_CODE: dict[str, tuple[PlanCache, CompileCache, CallLedger]] = {}


def _code_bundle(cache_dir: Path, ledger: CallLedger | None) -> tuple[PlanCache, CompileCache, CallLedger]:
    key = str(Path(cache_dir).resolve())
    hit = _CODE.get(key)
    if hit is not None:
        return hit
    led = ledger if ledger is not None else CallLedger()
    plans = PlanCache(led)
    cache = CompileCache(Path(key), led)
    _CODE[key] = (plans, cache, led)
    return plans, cache, led


class PhysicsSession:
    """One built DAG, one fused kernel, reusable contexts.

    Default _ctx is not concurrent-safe. Workers must spawn_context().
    """

    def __init__(
        self,
        spec: PlanSpec,
        plan: BuiltPlan,
        kernel: FusedKernel,
        ledger: CallLedger,
        constructor: str,
    ) -> None:
        self.spec = spec
        self.plan = plan
        self.kernel = kernel
        self.ledger = ledger
        self.constructor = constructor
        self._ctx = kernel.spawn()

    def spawn_context(self) -> FusedContext:
        return self.kernel.spawn()

    def source(self) -> str:
        return self.kernel.source

    def leaf_fortran(self) -> str:
        return self.kernel.leaf_fortran

    def _require_ctx(self, ctx: FusedContext | None) -> FusedContext:
        use = ctx or self._ctx
        if use.kernel_id != id(self.kernel):
            raise GuardError("G07", "context belongs to a different compiled kernel")
        if use.spec_id != self.spec.spec_id:
            raise GuardError("G07", "context spec_id mismatch")
        return use

    def evaluate(self, binding: RuntimeBinding, ctx: FusedContext | None = None) -> complex:
        value, _, _, _ = self.evaluate_timed(binding, ctx)
        return value

    def evaluate_timed(
        self, binding: RuntimeBinding, ctx: FusedContext | None = None
    ) -> tuple[complex, float, float, float]:
        """Returns value, validate_s, pack_s, kernel_s. Validate is billed, not omitted."""
        t0 = time.perf_counter()
        validate_binding(self.spec, binding)
        validate_s = time.perf_counter() - t0
        use = self._require_ctx(ctx)
        t1 = time.perf_counter()
        use.fill_sample(binding)
        pack_s = time.perf_counter() - t1
        t2 = time.perf_counter()
        out = use.run_filled(binding)
        kernel_s = time.perf_counter() - t2
        if self.ledger is not None:
            self.ledger.eval_calls[self.spec.spec_id] += 1
        return complex(out[0]), validate_s, pack_s, kernel_s

    def evaluate_batch(self, bindings: list[RuntimeBinding]) -> tuple[np.ndarray, list[str | None], float]:
        """Single-thread Fortran batch. Illegal samples rejected before kernel; not silent zeros."""
        n = len(bindings)
        values = np.full(n, np.nan + 0j, dtype=np.complex128)
        codes: list[str | None] = [None] * n
        legal: list[int] = []
        t0 = time.perf_counter()
        for i, b in enumerate(bindings):
            try:
                validate_binding(self.spec, b)
                legal.append(i)
            except GuardError as e:
                codes[i] = e.code
                if self.ledger is not None:
                    self.ledger.reject_calls[self.spec.spec_id] += 1
        if not legal:
            return values, codes, time.perf_counter() - t0
        ns = len(legal)
        nq = self.spec.n
        ntau = 2 * nq
        k = np.empty(ns, dtype=np.float64)
        q = np.empty((nq, ns), dtype=np.float64, order="F")
        tau = np.empty((ntau, ns), dtype=np.float64, order="F")
        t = np.empty(ns, dtype=np.float64)
        omega = np.empty(ns, dtype=np.float64)
        g = np.empty(ns, dtype=np.float64)
        L = np.empty(ns, dtype=np.float64)
        for j, i in enumerate(legal):
            b = bindings[i]
            k[j] = b.k
            q[:, j] = b.q
            tau[:, j] = b.tau
            t[j] = b.t
            omega[j] = b.omega
            g[j] = b.g
            L[j] = float(b.L)
        work = np.empty(max(1, self.kernel.ir.n_work), dtype=np.complex128)
        y = np.empty(ns, dtype=np.complex128)
        st = np.zeros(ns, dtype=np.int32)
        self.kernel.run_batch(k, q, tau, t, omega, g, L, work, y, st)
        for j, i in enumerate(legal):
            if st[j] != 0:
                codes[i] = f"FSTATUS_{int(st[j])}"
                values[i] = y[j]
            else:
                values[i] = y[j]
                if self.ledger is not None:
                    self.ledger.eval_calls[self.spec.spec_id] += 1
        return values, codes, time.perf_counter() - t0


def create_session(
    spec: PlanSpec,
    *,
    constructor: str = "R1",
    cache_dir: str | Path,
    ledger: CallLedger | None = None,
) -> PhysicsSession:
    """Build the spec once, compile fused leaf+DAG once, return a reusable session.

    constructor: "R1" or "CSE". CSE uses the same leaf Fortran writer as R1.
    Do not call this inside a sample loop. Code is cached; the session object is not shared.
    """
    if constructor not in _BUILDERS:
        raise KeyError(constructor)
    if spec.bands != 1 or spec.electron_interface != "scalar":
        raise ValueError("create_session fused kernel is scalar Holstein only")
    plans, cache, bundled_ledger = _code_bundle(Path(cache_dir), ledger)
    led = ledger if ledger is not None else bundled_ledger
    plan = plans.get(spec, builder=_BUILDERS[constructor])
    kernel = cache.compile_fused(spec, plan.graph, plan.leaves, tag=plan.report.constructor)
    return PhysicsSession(spec, plan, kernel, led, constructor)
