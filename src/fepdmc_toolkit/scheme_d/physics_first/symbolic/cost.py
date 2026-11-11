"""Call accounting. Unmeasured times are NOT_COMPUTED, never a silent 0.0."""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass, field

NOT_COMPUTED = "NOT_COMPUTED"


@dataclass
class CallLedger:
    """True invocation counts. Cache hits do not increment build/compile."""

    build_calls: dict[str, int] = field(default_factory=lambda: defaultdict(int))
    compile_calls: dict[str, int] = field(default_factory=lambda: defaultdict(int))
    reference_bind_calls: dict[str, int] = field(default_factory=lambda: defaultdict(int))
    packed_bind_calls: dict[str, int] = field(default_factory=lambda: defaultdict(int))
    eval_calls: dict[str, int] = field(default_factory=lambda: defaultdict(int))
    reject_calls: dict[str, int] = field(default_factory=lambda: defaultdict(int))

    def reset(self) -> None:
        for d in (
            self.build_calls,
            self.compile_calls,
            self.reference_bind_calls,
            self.packed_bind_calls,
            self.eval_calls,
            self.reject_calls,
        ):
            d.clear()

    def as_rows(self) -> list[dict]:
        keys = sorted(
            set(self.build_calls)
            | set(self.compile_calls)
            | set(self.eval_calls)
            | set(self.reject_calls)
            | set(self.reference_bind_calls)
            | set(self.packed_bind_calls)
        )
        rows = []
        for k in keys:
            rows.append(
                {
                    "spec_id": k,
                    "build_calls": self.build_calls.get(k, 0),
                    "compile_calls": self.compile_calls.get(k, 0),
                    "reference_bind_calls": self.reference_bind_calls.get(k, 0),
                    "packed_bind_calls": self.packed_bind_calls.get(k, 0),
                    "valid_evaluations": self.eval_calls.get(k, 0),
                    "rejects_before_kernel": self.reject_calls.get(k, 0),
                }
            )
        return rows


@dataclass
class BuildReport:
    spec_id: str
    constructor: str
    pairings_materialized: int
    recursive_calls: int
    unique_states: int
    memo_hits: int
    dag_nodes: int
    dag_edges: int
    n_leaves: int
    open_boundary_width: int
    build_s: float
    compile_s: str | float = NOT_COMPUTED
    bind_s: str | float = NOT_COMPUTED
    eval_s: str | float = NOT_COMPUTED
    peak_mb: str | float = NOT_COMPUTED
    notes: tuple[str, ...] = ()

    def as_row(self) -> dict:
        d = self.__dict__.copy()
        d["notes"] = "|".join(self.notes)
        return d
