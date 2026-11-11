"""Cost meter for recurrences. Hidden pairing materialization is a method fail."""

from __future__ import annotations

from dataclasses import dataclass, field


@dataclass
class CostMeter:
    constructor: str
    recursive_calls: int = 0
    unique_states: int = 0
    memo_hits: int = 0
    pairings_materialized: int = 0
    dag_nodes: int = 0
    dag_edges: int = 0
    leaf_bindings: int = 0
    build_s: float = 0.0
    compile_s: float = 0.0
    bind_s: float = 0.0
    eval_s: float = 0.0
    peak_mb: float = 0.0
    notes: tuple[str, ...] = field(default_factory=tuple)

    def as_row(self) -> dict:
        d = self.__dict__.copy()
        d["notes"] = "|".join(self.notes)
        return d
