"""Forward/reverse P1 score: occupancy reverse, not double-flip of Δt and ΔE."""

from __future__ import annotations

from fepdmc_toolkit.p1_da import logp_pkchange, reverse_logp


def test_legacy_double_flip_is_wrong():
    old, new, t_l, t_r = 0.0, 1.0, 0.0, 1.0
    forward = -(new - old) * (t_r - t_l)
    legacy = -(old - new) * (t_l - t_r)
    assert forward == -1.0
    assert legacy == -1.0
    assert legacy != -forward


def test_occupancy_reverse_is_antisymmetric():
    row = dict(tauL=1.0, tauR=2.0, w1=0.07, w2=0.12, tpl=-3.0, tpr=4.0)
    f = logp_pkchange(row, 0.3, 0.8)
    r = reverse_logp(row, 0.3, 0.8)
    assert abs(f + r) < 1e-14
