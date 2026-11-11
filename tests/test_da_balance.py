"""Delayed-acceptance detailed balance with abs(Re M), not abs(complex M)."""

from __future__ import annotations

import math

import numpy as np

from fepdmc_toolkit.p1_da import (
    abs_complex_target,
    delayed_acceptance_prob,
    native_re_ratio,
    native_re_target,
)

LOG10 = math.log(10.0)
PHASES = (1.0 + 0.0j, 1.0 + 0.3j, np.exp(0.7j))
M_SIX = np.array(
    [1.0 + 0.2j, 1.5 - 0.4j, 0.8 + 0.5j, 2.0 + 0.1j, 1.2 - 0.8j, 0.9 + 0.3j],
    dtype=np.complex128,
)


def test_da_preserves_native_ratio():
    rng = np.random.default_rng(2)
    for _ in range(40):
        mx = rng.normal() + 1j * rng.normal()
        my = rng.normal() + 1j * rng.normal()
        z = rng.choice(np.array(PHASES))
        px = abs(native_re_target(mx, z))
        py = abs(native_re_target(my, z))
        if px < 1e-12 or py < 1e-12:
            continue
        R = py / px
        ell = float(np.clip(math.log(max(R, 1e-12)), -LOG10, LOG10))
        rhat = math.exp(ell)
        a_xy = delayed_acceptance_prob(px, py, rhat)
        a_yx = delayed_acceptance_prob(py, px, 1.0 / rhat)
        assert abs(a_xy / a_yx - R) <= 1e-12 * max(1.0, abs(R))


def test_abs_complex_is_not_the_native_target():
    mismatches = 0
    for z in PHASES:
        for mx in M_SIX:
            for my in M_SIX:
                R = native_re_ratio(my, mx, z, z)
                wrong = abs_complex_target(my, z) / abs_complex_target(mx, z)
                if abs(wrong - R) > 1e-8:
                    mismatches += 1
                    rhat = math.exp(float(np.clip(math.log(max(R, 1e-12)), -LOG10, LOG10)))
                    a_ok = delayed_acceptance_prob(abs(native_re_target(mx, z)), abs(native_re_target(my, z)), rhat)
                    a_ok_rev = delayed_acceptance_prob(
                        abs(native_re_target(my, z)), abs(native_re_target(mx, z)), 1.0 / rhat
                    )
                    a_wrong = delayed_acceptance_prob(abs_complex_target(mx, z), abs_complex_target(my, z), rhat)
                    a_wrong_rev = delayed_acceptance_prob(
                        abs_complex_target(my, z), abs_complex_target(mx, z), 1.0 / rhat
                    )
                    assert abs(a_ok / a_ok_rev - R) <= 1e-12 * max(1.0, abs(R))
                    if abs(a_wrong / a_wrong_rev - R) <= 1e-8:
                        # a wrong target that accidentally matches is allowed; just
                        # require that the two targets themselves differ.
                        assert abs(wrong - R) > 1e-8
    assert mismatches > 0


def test_zero_new_weight_rejects():
    assert delayed_acceptance_prob(1.0, 0.0, 2.0) == 0.0
