"""S05: Python leaf evaluator. Must not call builder or compiler."""

from __future__ import annotations

import time

import numpy as np

from .cost import CallLedger
from .guards import GuardError, RuntimeBinding, validate_binding
from .leaves import LeafTable, eval_leaf
from .plan_spec import PlanSpec


class ReferenceBinder:
    def __init__(self, spec: PlanSpec, leaves: LeafTable, ledger: CallLedger | None = None) -> None:
        self.spec = spec
        self.leaves = leaves
        self.ledger = ledger

    def bind(self, binding: RuntimeBinding) -> dict[str, np.ndarray]:
        validate_binding(self.spec, binding)
        if self.ledger is not None:
            self.ledger.reference_bind_calls[self.spec.spec_id] += 1
        env: dict[str, np.ndarray] = {}
        for name in self.leaves.order:
            env[name] = eval_leaf(self.leaves.descriptors[name], binding)
        return env

    def bind_timed(self, binding: RuntimeBinding) -> tuple[dict[str, np.ndarray], float]:
        t0 = time.perf_counter()
        env = self.bind(binding)
        return env, time.perf_counter() - t0


def try_bind(binder: ReferenceBinder, binding: RuntimeBinding) -> tuple[dict[str, np.ndarray] | None, str | None]:
    try:
        return binder.bind(binding), None
    except GuardError as e:
        if binder.ledger is not None:
            binder.ledger.reject_calls[binder.spec.spec_id] += 1
        return None, e.code
