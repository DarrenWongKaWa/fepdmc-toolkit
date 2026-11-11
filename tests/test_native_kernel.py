"""Stage 4: NativeKernelIR lowering, interpreter, laziness, mutations."""

from __future__ import annotations

import math

import numpy as np
import pytest

from fepdmc_toolkit.diagram_compiler import (
    Binding,
    build_scalar_ir,
    build_twoband_ir,
    compile_delayed_acceptance,
    interpret_native_kernel,
    lower_native_kernel,
    lower_native_score_kernel,
    native_acceptance_probability,
    validate_native_ir,
)
from fepdmc_toolkit.diagram_compiler.evaluator import torus_point
from fepdmc_toolkit.diagram_compiler.native_ir import NativeIRError as IRErr
from fepdmc_toolkit.diagram_compiler.native_validate import exact_candidate_ops_on_reject_path


def _binding(n: int, *, model="scalar", variant="x") -> Binding:
    L = 4
    if variant == "x":
        k_idx, q0, tau_step, t, omega, g = 0, 1, 0.05, 1.0, 0.8, 0.5
    else:
        k_idx, q0, tau_step, t, omega, g = 1, 2, 0.07, 1.1, 0.9, 0.4
    return Binding(
        k=torus_point(k_idx, L),
        q=tuple(torus_point(q0 + i, L) for i in range(n)),
        tau=tuple(tau_step * (i + 1) for i in range(2 * n)),
        t=t,
        omega=omega,
        g=g,
        L=L,
        delta=0.35,
        gap=0.5,
    )


def test_binding_layout_two_band_n1():
    from fepdmc_toolkit.diagram_compiler.native_ir import binding_layout_for

    layout = binding_layout_for("Holstein_twoband_hermitian", 1)
    names = [f.name for f in layout.fields]
    assert names[:7] == ["k", "q", "tau", "t", "omega", "g", "L"]
    assert "delta" in names and "gap" in names
    q = next(f for f in layout.fields if f.name == "q")
    tau = next(f for f in layout.fields if f.name == "tau")
    assert q.shape == (1,) and tau.shape == (2,)


def test_twoband_kernel_round_trip_and_validate():
    kernel = compile_delayed_acceptance(build_twoband_ir(1))
    nir = lower_native_kernel(kernel)
    validate_native_ir(nir)
    revived = type(nir).from_json(nir.to_json())
    assert revived.to_dict() == nir.to_dict()
    assert exact_candidate_ops_on_reject_path(nir) == []


def test_exact_y_absent_from_stage1_reject_path():
    nir = lower_native_kernel(compile_delayed_acceptance(build_twoband_ir(1)))
    assert exact_candidate_ops_on_reject_path(nir) == []
    reject_ops = [op.op for block in nir.blocks if block.id in ("entry", "reject_stage1") for op in block.ops]
    assert "call_graph" in reject_ops
    assert all(
        op.attrs.get("graph") != "exact"
        for block in nir.blocks
        if block.id in ("entry", "reject_stage1")
        for op in block.ops
    )


def _compare(kernel, nir, x, y, log_q, u1, u2, **cache):
    py = kernel.evaluate_transition(x, y, log_q, u1, u2, **cache)
    native = interpret_native_kernel(
        nir,
        x=x,
        y=y,
        log_q=log_q,
        u1=u1,
        u2=u2,
        exact_x_valid=cache.get("exact_x") is not None,
        exact_log_weight_x=cache.get("exact_log_weight_x", 0.0),
    )
    assert native.status == 0
    assert native.accepted == py.accepted
    assert native.exact_y_evaluated == py.exact_y_evaluated
    assert native.stage1_pass == py.stage1_pass
    assert native.stage2_pass == py.stage2_pass
    assert native.ell_hat == pytest.approx(py.ell_hat)
    if py.ell_R is None:
        assert native.ell_R is None
    else:
        assert native.ell_R == pytest.approx(py.ell_R)


@pytest.mark.parametrize("n", [1, 2, 3, 4])
def test_scalar_differential_three_paths(n):
    kernel = compile_delayed_acceptance(build_scalar_ir(n))
    nir = lower_native_kernel(kernel)
    x, y = _binding(n, variant="x"), _binding(n, variant="y")
    probe = kernel.evaluate_transition(x, y, 0.0, 1e-16, 1e-16)
    alpha2 = 1.0 if probe.ell_R - probe.ell_hat >= 0 else math.exp(probe.ell_R - probe.ell_hat)
    u2_reject = min(0.999999, (alpha2 + 1.0) / 2.0)
    _compare(kernel, nir, x, y, 0.0, 0.9, 0.1)
    if alpha2 < 1.0:
        _compare(kernel, nir, x, y, 0.0, 1e-16, u2_reject)
    _compare(kernel, nir, x, y, 0.0, 1e-16, 1e-16)
    _compare(kernel, nir, x, y, 0.4, 1e-16, 1e-16)


def test_twoband_n1_differential_and_laziness():
    kernel = compile_delayed_acceptance(build_twoband_ir(1))
    nir = lower_native_kernel(kernel)
    x, y = _binding(1, variant="x"), _binding(1, variant="y")
    _compare(kernel, nir, x, y, 0.0, 0.9, 0.1)
    native_reject = interpret_native_kernel(nir, x=x, y=y, log_q=0.0, u1=0.9, u2=0.1)
    assert native_reject.exact_y_evaluated is False
    probe = kernel.evaluate_transition(x, y, 0.0, 1e-16, 1e-16)
    alpha2 = math.exp(min(0.0, probe.ell_R - probe.ell_hat))
    _compare(kernel, nir, x, y, 0.0, 1e-16, (alpha2 + 1.0) / 2.0)
    _compare(kernel, nir, x, y, 0.0, 1e-16, 1e-16)


def test_cached_exact_x_matches_python():
    kernel = compile_delayed_acceptance(build_twoband_ir(1))
    nir = lower_native_kernel(kernel)
    x, y = _binding(1, variant="x"), _binding(1, variant="y")
    log_px = kernel.exact_log_weight(x)
    _compare(kernel, nir, x, y, 0.0, 1e-16, 1e-16, exact_log_weight_x=log_px, exact_x=x)


def test_finite_state_via_native_scores():
    pi = np.array([0.5, 0.3, 0.2])
    cheap = np.array([0.1, 0.55, 0.35])
    q = np.array([[0.0, 0.8, 0.2], [0.3, 0.0, 0.7], [0.6, 0.4, 0.0]])
    n = 3
    p = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            log_q = math.log(q[j, i]) - math.log(q[i, j])
            a = native_acceptance_probability(
                math.log(cheap[i]), math.log(cheap[j]), math.log(pi[i]), math.log(pi[j]), log_q
            )
            # Cross-check the score kernel path with always-pass uniforms.
            score = lower_native_score_kernel()
            result = interpret_native_kernel(
                score,
                log_q=log_q,
                u1=1e-16,
                u2=1e-16,
                log_wx=math.log(cheap[i]),
                log_wy=math.log(cheap[j]),
                log_px=math.log(pi[i]),
                log_py=math.log(pi[j]),
            )
            assert result.ell_hat == pytest.approx(math.log(cheap[j]) - math.log(cheap[i]))
            p[i, j] = q[i, j] * a
        p[i, i] = 1.0 - p[i].sum()
    residual = max(abs(pi[i] * p[i, j] - pi[j] * p[j, i]) for i in range(n) for j in range(n))
    assert residual <= 1e-12
    assert np.allclose(pi @ p, pi, atol=1e-12, rtol=1e-10)


def test_mutation_exact_on_entry_fails_validator():
    nir = lower_native_kernel(compile_delayed_acceptance(build_twoband_ir(1)))
    entry = nir.block("entry")
    from fepdmc_toolkit.diagram_compiler.native_ir import REGION_EXACT, TYPE_C128, NativeBlock, NativeOp

    bad_ops = entry.ops + (
        NativeOp("sneak", "call_graph", TYPE_C128, (), {"graph": "exact", "binding": "y"}, REGION_EXACT),
    )
    blocks = []
    for block in nir.blocks:
        if block.id == "entry":
            blocks.append(NativeBlock("entry", bad_ops, entry.term))
        else:
            blocks.append(block)
    import dataclasses

    mutated = dataclasses.replace(nir, blocks=tuple(blocks))
    with pytest.raises(IRErr):
        validate_native_ir(mutated)


def test_undefined_id_and_bad_branch_fail():
    import dataclasses

    from fepdmc_toolkit.diagram_compiler.native_ir import NativeBlock, NativeTerminator

    nir = lower_native_score_kernel()
    bad = dataclasses.replace(
        nir,
        blocks=tuple(
            NativeBlock(b.id, b.ops, NativeTerminator("br", (), {"target": "missing"})) if b.id == "accept" else b
            for b in nir.blocks
        ),
    )
    with pytest.raises(IRErr):
        validate_native_ir(bad)


def test_interpreter_does_not_call_python_da_kernel(monkeypatch):
    import fepdmc_toolkit.diagram_compiler.da_kernel as da

    def boom(*_a, **_k):
        raise AssertionError("evaluate_transition used as shortcut")

    monkeypatch.setattr(da.DelayedAcceptanceKernel, "evaluate_transition", boom)
    kernel = compile_delayed_acceptance(build_twoband_ir(1))
    nir = lower_native_kernel(kernel)
    x, y = _binding(1, variant="x"), _binding(1, variant="y")
    result = interpret_native_kernel(nir, x=x, y=y, log_q=0.0, u1=1e-16, u2=1e-16)
    assert result.accepted is True
