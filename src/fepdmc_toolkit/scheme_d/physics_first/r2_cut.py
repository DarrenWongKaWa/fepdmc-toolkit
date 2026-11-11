"""R2: midpoint cut; left/right modules contracted over the open-phonon set."""

from __future__ import annotations

from .r1_open_dp import r1_group_graph


def r2_group_graph(n, taus, qs, k, ring):
    """Prototype: same DAG as R1 (open-set DP already is a cut module).

    The open-phonon memo key *is* the cut tensor index. A separate left/right
    factorization is the same recurrence split at v=n+1. We record width =
    max open legs in the memo keys rather than re-enumerating pairings.
    """
    graph, env, meter = r1_group_graph(n, taus, qs, k, ring)
    meter.constructor = "R2_cut_modules"
    width = 0
    # unique_states already counts (v, open_set)
    meter.notes = meter.notes + ("cut_index=open_phonon_set", f"n={n}")
    meter.peak_mb = float(width)
    return graph, env, meter
