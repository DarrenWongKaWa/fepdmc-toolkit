"""M3: pairwise delayed-acceptance algebra for Design B.

A_xy / A_yx = R_exact when ell_hat is reciprocal and ell_R includes Hastings.
"""

from __future__ import annotations

import math

import pytest

from fepdmc_toolkit.diagram_compiler.da_kernel import (
    acceptance_probability,
    delayed_acceptance_identity,
    ell_hat_state_weight,
    ell_R_design_b,
)


def test_ell_hat_is_antisymmetric():
    log_wx, log_wy = math.log(0.2), math.log(0.8)
    fwd = ell_hat_state_weight(log_wx, log_wy)
    rev = ell_hat_state_weight(log_wy, log_wx)
    assert rev == pytest.approx(-fwd)
    assert fwd == pytest.approx(math.log(0.8) - math.log(0.2))


def test_design_b_ell_R_includes_hastings():
    log_px, log_py = math.log(0.5), math.log(0.1)
    log_q_rev_minus_fwd = math.log(0.3) - math.log(0.9)
    ell = ell_R_design_b(log_px, log_py, log_q_rev_minus_fwd)
    assert ell == pytest.approx((log_py - log_px) + log_q_rev_minus_fwd)


@pytest.mark.parametrize(
    "wx,wy,px,py,qxy,qyx",
    [
        (0.4, 0.2, 0.5, 0.5, 0.5, 0.5),
        (0.1, 0.9, 0.8, 0.2, 0.5, 0.5),
        (0.3, 0.3, 0.1, 0.9, 0.2, 0.8),
        (0.05, 0.7, 0.4, 0.6, 0.9, 0.1),
        (1.5, 0.2, 0.25, 0.75, 0.4, 0.6),
    ],
)
def test_design_b_identity_holds(wx, wy, px, py, qxy, qyx):
    ell_hat = ell_hat_state_weight(math.log(wx), math.log(wy))
    log_q = math.log(qyx) - math.log(qxy)
    ell_R = ell_R_design_b(math.log(px), math.log(py), log_q)
    a_xy = acceptance_probability(ell_hat, ell_R)
    a_yx = acceptance_probability(-ell_hat, -ell_R)
    r_exact = (py * qyx) / (px * qxy)
    assert delayed_acceptance_identity(a_xy, a_yx, r_exact)


def test_removing_stage2_breaks_identity():
    ell_hat = ell_hat_state_weight(math.log(0.1), math.log(0.9))
    a_xy = math.exp(min(0.0, ell_hat))  # stage1 only
    a_yx = math.exp(min(0.0, -ell_hat))
    r_exact = 0.2 / 0.8
    assert not delayed_acceptance_identity(a_xy, a_yx, r_exact)


def test_dropping_minus_ell_hat_from_stage2_breaks_identity():
    ell_hat = ell_hat_state_weight(math.log(0.1), math.log(0.9))
    ell_R = ell_R_design_b(math.log(0.8), math.log(0.2), 0.0)
    a_xy = min(1.0, math.exp(ell_hat)) * min(1.0, math.exp(ell_R))
    a_yx = min(1.0, math.exp(-ell_hat)) * min(1.0, math.exp(-ell_R))
    r_exact = 0.2 / 0.8
    assert not delayed_acceptance_identity(a_xy, a_yx, r_exact)
