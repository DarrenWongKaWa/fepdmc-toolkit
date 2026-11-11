"""R3: legal dummy-variable canonicalization vs illegal bound-q merge."""

from __future__ import annotations

from .r1_open_dp import r1_group_graph


def r3_legal_bound_x(n, taus, qs, k, ring):
    """Bound x: phonon ids are not dummies. Legal key keeps identities (R1)."""
    graph, env, meter = r1_group_graph(n, taus, qs, k, ring)
    meter.constructor = "R3_legal_bound_ids"
    meter.notes = meter.notes + ("bound_q_not_renamed",)
    return graph, env, meter


def illegal_merge_key_open_count_only(n: int, open_ph, v: int) -> tuple:
    """DIAGNOSTIC_WRONG_MODEL: drops q identity. Must not be used for bound x."""
    return (v, len(open_ph))
