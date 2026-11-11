"""Exact/cheap/proposal independence of the DA kernel runtime."""

from __future__ import annotations

from fepdmc_toolkit.diagram_compiler import Binding, build_twoband_ir, compile_delayed_acceptance
from fepdmc_toolkit.diagram_compiler.evaluator import torus_point


def _binding(variant: str = "x") -> Binding:
    L = 4
    if variant == "x":
        return Binding(
            k=torus_point(0, L),
            q=(torus_point(1, L),),
            tau=(0.05, 0.12),
            t=1.0,
            omega=0.8,
            g=0.5,
            L=L,
            delta=0.35,
            gap=0.5,
        )
    return Binding(
        k=torus_point(1, L),
        q=(torus_point(2, L),),
        tau=(0.07, 0.19),
        t=1.1,
        omega=0.9,
        g=0.4,
        L=L,
        delta=0.35,
        gap=0.5,
    )


def test_kernel_does_not_import_r1_or_p1_at_runtime(monkeypatch):
    import fepdmc_toolkit.diagram_compiler.da_kernel as da
    import fepdmc_toolkit.diagram_compiler.evaluator as ev
    import fepdmc_toolkit.diagram_compiler.graph as graph

    ir = build_twoband_ir(1)
    kernel = compile_delayed_acceptance(ir)

    def boom(*_a, **_k):
        raise AssertionError("forbidden backend")

    monkeypatch.setattr(ev, "evaluate_diagramwise", boom)
    monkeypatch.setattr(graph, "lower_grouped", boom)
    # compile_evaluator / evaluate_dag are used by the exact wrapper created
    # at compile time; the kernel must not re-lower or call R1.
    monkeypatch.setattr(graph, "lower_diagramwise", boom)
    result = kernel.evaluate_transition(_binding("x"), _binding("y"), 0.0, 1e-16, 1e-16)
    assert result.ell_hat == result.ell_hat
    assert "scheme_d" not in da.__dict__
    assert "linear_da" not in da.__dict__


def test_proposal_does_not_recompute_diagram_weights():
    ir = build_twoband_ir(1)
    kernel = compile_delayed_acceptance(ir)
    x, y = _binding("x"), _binding("y")
    a = kernel.evaluate_transition(x, y, 0.0, 1e-16, 1e-16)
    b = kernel.evaluate_transition(x, y, 0.4, 1e-16, 1e-16)
    assert a.ell_hat == b.ell_hat
    assert a.cheap_log_weight_x == b.cheap_log_weight_x
    assert a.ell_R != b.ell_R
