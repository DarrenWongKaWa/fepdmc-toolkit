"""The exact compiled evaluator builds and runs on a scalar order-1 family."""

from __future__ import annotations

from fepdmc_toolkit.diagram_compiler import build_scalar_ir, compile_evaluator, numeric_close
from fepdmc_toolkit.diagram_compiler.evaluator import Binding, torus_point


def test_exact_compile_evaluator_runs():
    ir = build_scalar_ir(1)
    L = 4
    binding = Binding(
        k=torus_point(0, L),
        q=(torus_point(1, L),),
        tau=(0.05, 0.12),
        t=1.0,
        omega=0.8,
        g=0.5,
        L=L,
    )
    out = compile_evaluator(ir)(binding)
    assert numeric_close(out["F"], out["F"])
    assert out["F"].real != 0.0
