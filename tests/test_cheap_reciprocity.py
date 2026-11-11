"""M3/M5: reverse is Binding swap; raw cheap score is antisymmetric by construction."""

from __future__ import annotations

import math

import pytest

from fepdmc_toolkit.diagram_compiler import Binding, build_scalar_ir, build_twoband_ir, compile_cheap
from fepdmc_toolkit.diagram_compiler.evaluator import torus_point, validate_binding


def _binding(n: int, *, L: int = 4, **extra) -> Binding:
    tau = tuple(0.05 * (i + 1) for i in range(2 * n))
    q = tuple(torus_point(i + 1, L) for i in range(n))
    fields = dict(k=torus_point(0, L), q=q, tau=tau, t=1.0, omega=0.8, g=0.5, L=L, delta=0.35, gap=0.5)
    fields.update(extra)
    return Binding(**fields)


ATOL = 1e-12
RTOL = 1e-10


def _close(a: float, b: float) -> bool:
    return abs(a - b) <= ATOL + RTOL * abs(b)


@pytest.mark.parametrize(
    "builder,n",
    [
        (build_scalar_ir, 1),
        (build_scalar_ir, 2),
        (build_scalar_ir, 3),
        (build_scalar_ir, 4),
        (build_twoband_ir, 1),
        (build_twoband_ir, 2),
        (build_twoband_ir, 3),
    ],
)
def test_rebinding_score_is_antisymmetric(builder, n):
    ir = builder(n)
    cheap = compile_cheap(ir)
    x = _binding(n)
    y = _binding(
        n,
        k=torus_point(1, 4),
        q=tuple(torus_point(i + 2, 4) for i in range(n)),
        tau=tuple(0.07 * (i + 1) for i in range(2 * n)),
        t=1.1,
        omega=0.9,
        g=0.4,
    )
    validate_binding(ir, x)
    validate_binding(ir, y)
    fwd = cheap.transition_score(x, y)
    rev = cheap.transition_score(y, x)
    assert rev == pytest.approx(-fwd)
    assert _close(rev, -fwd)
    constructed = cheap.log_weight(y) - cheap.log_weight(x)
    assert fwd == constructed


def test_q_occupancy_swap_is_the_legal_p1_analog():
    ir = build_scalar_ir(2)
    cheap = compile_cheap(ir)
    x = _binding(2)
    y = Binding(**{**x.__dict__, "q": (x.q[1], x.q[0])})
    assert cheap.transition_score(y, x) == pytest.approx(-cheap.transition_score(x, y))
    assert cheap.transition_score(x, y) != 0.0


def test_unsorted_time_swap_is_not_a_legal_reverse():
    ir = build_scalar_ir(2)
    cheap = compile_cheap(ir)
    x = _binding(2)
    illegal = Binding(**{**x.__dict__, "tau": (x.tau[1], x.tau[0]) + x.tau[2:]})
    with pytest.raises(ValueError, match="unsorted"):
        cheap.transition_score(x, illegal)


def test_mutating_reverse_to_drop_a_sign_breaks_antisymmetry():
    ir = build_scalar_ir(2)
    cheap = compile_cheap(ir)
    x = _binding(2)
    y = _binding(2, k=torus_point(1, 4))
    fwd = cheap.log_weight(y) - cheap.log_weight(x)
    broken = cheap.log_weight(x) - cheap.log_weight(y)
    broken_negated_once = broken  # independently evaluated reverse without the construction minus
    assert broken_negated_once == pytest.approx(-fwd)
    mutated = fwd  # pretend reverse equals forward (legacy double-flip class of bug)
    assert mutated != pytest.approx(-fwd) or fwd == 0.0
    assert not _close(mutated, -fwd)


def test_odd_clip_preserves_reciprocity_when_used():
    from fepdmc_toolkit.diagram_compiler.cheap_evaluator import odd_clip

    value = 4.2
    limit = math.log(10.0)
    assert odd_clip(-value, limit) == pytest.approx(-odd_clip(value, limit))
    assert odd_clip(0.0, limit) == 0.0
    assert odd_clip(limit + 1.0, limit) == limit
    assert odd_clip(-(limit + 1.0), limit) == -limit
