"""F11: common recurrence interface. Unique writer for physics-first DAGs."""

from __future__ import annotations

from .baselines import enum_sum_graph
from .r1_open_dp import r1_group_graph
from .r2_cut import r2_group_graph
from .r3_canonical import r3_legal_bound_x
from .r4_graded import r4_nested_graph

# R1 here is the historical numeric-leaf builder (float k_eff/dtau names).
# Compile-once reuse lives in physics_first.symbolic, not this map.
CONSTRUCTORS = {
    "R1": r1_group_graph,
    "R2": r2_group_graph,
    "R3": r3_legal_bound_x,
    "R4": r4_nested_graph,
    "B-enum-CSE": lambda n, taus, qs, k, ring: enum_sum_graph(n, taus, qs, k, ring, hash_cons=True),
    "B-enum": lambda n, taus, qs, k, ring: enum_sum_graph(n, taus, qs, k, ring, hash_cons=False),
}


def build(name: str, n, taus, qs, k, ring):
    """Historical per-sample numeric-leaf builder. Not compile-once reuse."""
    if name not in CONSTRUCTORS:
        raise KeyError(name)
    return CONSTRUCTORS[name](n, taus, qs, k, ring)


def create_r1_session(spec, cache_dir, ledger=None):
    """Public compile-once session. Keep the object; call evaluate per sample."""
    from .symbolic.session import create_session

    return create_session(spec, constructor="R1", cache_dir=cache_dir, ledger=ledger)
