"""M8: P1 correspondence is a research comparison, not an identity."""

from __future__ import annotations

import math

from fepdmc_toolkit.diagram_compiler import Binding, build_scalar_ir, compile_cheap
from fepdmc_toolkit.diagram_compiler.evaluator import torus_point


def _xi(k: float, t: float) -> float:
    return 2.0 * t * (1.0 - math.cos(k))


def test_p1_local_gel_ratio_is_not_the_grouped_state_weight_ratio():
    """Historical P1 ell_hat = log P_kchange on one occupancy swap.

    Cheap stage-2 score is log W_hat(y) - log W_hat(x) of the grouped F_n
    propagator projection. These are different objects.
    """
    ir = build_scalar_ir(2)
    cheap = compile_cheap(ir)
    L = 4
    x = Binding(
        k=torus_point(0, L),
        q=(torus_point(1, L), torus_point(2, L)),
        tau=(0.05, 0.12, 0.20, 0.31),
        t=1.0,
        omega=0.8,
        g=0.5,
        L=L,
    )
    y = Binding(**{**x.__dict__, "q": (x.q[1], x.q[0])})
    grouped = cheap.transition_score(x, y)
    # P1-like local mid-segment Gel ratio for interval (tau0, tau1) using
    # k_eff = k - q_open. This is a hand local ratio, not the compiler.
    dtau = x.tau[1] - x.tau[0]
    gold = math.exp(-_xi(x.k - x.q[0], x.t) * dtau)
    gnew = math.exp(-_xi(y.k - y.q[0], y.t) * dtau)
    local = math.log(gnew / gold)
    assert grouped != 0.0
    assert abs(grouped - local) > 1e-8


def test_p1_correspondence_classification_is_explicit():
    from fepdmc_toolkit.diagram_compiler.cheap_evaluator import P1_CORRESPONDENCE

    assert P1_CORRESPONDENCE["verdict"] == "NOT_IDENTICAL"
    assert "representation_mismatch" in P1_CORRESPONDENCE["classes"]
    assert "proposal_metadata_missing" in P1_CORRESPONDENCE["classes"]
    assert "compiler_policy_different" in P1_CORRESPONDENCE["classes"]
    assert P1_CORRESPONDENCE["cheap_semantics"] == "STATE_WEIGHT"
    assert P1_CORRESPONDENCE["p1_semantics"] == "TRANSITION_SCORE"
