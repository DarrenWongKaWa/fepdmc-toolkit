"""Stage 5: NativeKernelIR Fortran backend vs independent interpreter."""

from __future__ import annotations

import ast
import hashlib
import math
from pathlib import Path

import numpy as np
import pytest

from fepdmc_toolkit.diagram_compiler import (
    Binding,
    NativeIRError,
    build_scalar_ir,
    build_twoband_ir,
    compile_delayed_acceptance,
    compile_fortran_source,
    emit_driver,
    emit_fortran,
    interpret_native_kernel,
    lower_native_kernel,
    lower_native_score_kernel,
)
from fepdmc_toolkit.diagram_compiler.evaluator import torus_point
from fepdmc_toolkit.diagram_compiler.fortran_codegen import BACKEND, native_ir_digest, validate_generated_fortran
from fepdmc_toolkit.diagram_compiler.fortran_compile import ATOL, RTOL, compile_native_kernel, gfortran_available
from fepdmc_toolkit.diagram_compiler.native_ir import TYPE_C128, NativeOp

ROOT = Path(__file__).resolve().parents[1]
COMPILER = ROOT / "src/fepdmc_toolkit/diagram_compiler"


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


def _match(py, ft):
    assert ft.status == py.status
    if py.status != 0:
        return
    assert ft.accepted == py.accepted
    assert ft.reject_stage == py.reject_stage
    assert ft.exact_y_evaluated == py.exact_y_evaluated
    assert ft.ell_hat == pytest.approx(py.ell_hat, abs=ATOL, rel=RTOL)
    if py.ell_R is None:
        assert ft.ell_R_valid is False
    else:
        assert ft.ell_R_valid is True
        assert ft.ell_R == pytest.approx(py.ell_R, abs=ATOL, rel=RTOL)


def _paths(kernel, x, y):
    probe = kernel.evaluate_transition(x, y, 0.0, 1e-16, 1e-16)
    alpha2 = 1.0 if probe.ell_R - probe.ell_hat >= 0 else math.exp(probe.ell_R - probe.ell_hat)
    u2_reject = min(0.999999, (alpha2 + 1.0) / 2.0)
    cases = [("stage1_reject", 0.0, 0.9, 0.1)]
    if probe.ell_hat < 0:
        pass
    else:
        cases = [("stage1_reject", 0.0, 0.9, 0.1)]
    if alpha2 < 1.0:
        cases.append(("stage2_reject", 0.0, 1e-16, u2_reject))
    cases.append(("accept", 0.0, 1e-16, 1e-16))
    cases.append(("asymmetric_q", 0.4, 1e-16, 1e-16))
    return cases


pytestmark = pytest.mark.skipif(not gfortran_available(), reason="gfortran not found")


@pytest.fixture(scope="module")
def twoband():
    kernel = compile_delayed_acceptance(build_twoband_ir(1))
    nir = lower_native_kernel(kernel)
    art = compile_native_kernel(nir)
    return kernel, nir, art


@pytest.fixture(scope="module")
def score_art():
    nir = lower_native_score_kernel()
    return nir, compile_native_kernel(nir)


def test_backend_imports_native_ir_only():
    tree = ast.parse((COMPILER / "fortran_codegen.py").read_text())
    imported = []
    for node in tree.body:
        if isinstance(node, ast.ImportFrom):
            imported.append(node.module)
        elif isinstance(node, ast.Import):
            imported.extend(alias.name for alias in node.names)
    joined = " ".join(imported or [])
    for banned in (
        "ir",
        "cheap_policy",
        "cheap_evaluator",
        "da_kernel",
        "canonical",
        "evaluator",
        "proposal",
        "target_policy",
        "r1_import",
    ):
        assert banned not in imported
        assert f".{banned}" not in joined
    assert any(name and name.endswith("native_ir") for name in imported)


def test_emit_refuses_unknown_backend_and_opcode(twoband):
    _kernel, nir, _art = twoband
    with pytest.raises(NativeIRError):
        emit_fortran(nir, backend="c_v1")
    import dataclasses

    bad_ops = nir.block("entry").ops + (NativeOp("bogus", "prim.not_a_thing", TYPE_C128),)
    from fepdmc_toolkit.diagram_compiler.native_ir import NativeBlock

    blocks = tuple(NativeBlock(b.id, bad_ops, b.term) if b.id == "entry" else b for b in nir.blocks)
    mutated = dataclasses.replace(nir, blocks=blocks)
    with pytest.raises(NativeIRError):
        emit_fortran(mutated)


def test_source_determinism(twoband):
    _kernel, nir, _art = twoband
    a = emit_fortran(nir)
    b = emit_fortran(nir)
    assert a == b
    assert hashlib.sha256(a.encode()).hexdigest() == hashlib.sha256(b.encode()).hexdigest()
    assert "20" not in a.split("native_ir_digest=")[0] or "timestamp" not in a.lower()
    assert "datetime" not in a.lower()


def test_twoband_n1_three_paths_and_laziness(twoband):
    kernel, nir, art = twoband
    x, y = _binding(1, variant="x"), _binding(1, variant="y")
    for name, log_q, u1, u2 in _paths(kernel, x, y):
        py = interpret_native_kernel(nir, x=x, y=y, log_q=log_q, u1=u1, u2=u2)
        ft = art.run_diagram(x, y, log_q, u1, u2, layout=nir.layout)
        _match(py, ft)
        if name == "stage1_reject":
            assert py.exact_y_evaluated is False
            assert ft.n_exact_graph == 0
            assert ft.n_electron_twoband == 0
            assert ft.n_vertex_sigmaz == 0
            assert ft.n_matmul == 0
            assert ft.n_trace == 0
        else:
            assert ft.n_exact_graph > 0
            assert ft.n_electron_twoband > 0
            assert ft.n_trace > 0


def test_cached_exact_x_skips_exact_x_graph(twoband):
    kernel, nir, art = twoband
    x, y = _binding(1, variant="x"), _binding(1, variant="y")
    log_px = kernel.exact_log_weight(x)
    py = interpret_native_kernel(
        nir, x=x, y=y, log_q=0.0, u1=1e-16, u2=1e-16, exact_x_valid=True, exact_log_weight_x=log_px
    )
    ft = art.run_diagram(x, y, 0.0, 1e-16, 1e-16, exact_x_valid=True, exact_log_weight_x=log_px, layout=nir.layout)
    _match(py, ft)
    assert ft.n_exact_graph == 1
    uncached = art.run_diagram(x, y, 0.0, 1e-16, 1e-16, layout=nir.layout)
    assert uncached.n_exact_graph == 2


@pytest.mark.parametrize("n", [1, 2, 3, 4])
def test_scalar_orders(n, tmp_path):
    kernel = compile_delayed_acceptance(build_scalar_ir(n))
    nir = lower_native_kernel(kernel)
    art = compile_native_kernel(nir, workdir=tmp_path / f"s{n}")
    x, y = _binding(n, variant="x"), _binding(n, variant="y")
    for name, log_q, u1, u2 in _paths(kernel, x, y):
        py = interpret_native_kernel(nir, x=x, y=y, log_q=log_q, u1=u1, u2=u2)
        ft = art.run_diagram(x, y, log_q, u1, u2, layout=nir.layout)
        _match(py, ft)
        if name == "stage1_reject":
            assert ft.n_exact_graph == 0
        else:
            assert ft.n_exact_graph > 0


def test_twoband_n2_if_positive(tmp_path):
    kernel = compile_delayed_acceptance(build_twoband_ir(2))
    nir = lower_native_kernel(kernel)
    art = compile_native_kernel(nir, workdir=tmp_path / "tb2")
    x, y = _binding(2, variant="x"), _binding(2, variant="y")
    py = interpret_native_kernel(nir, x=x, y=y, log_q=0.0, u1=1e-16, u2=1e-16)
    ft = art.run_diagram(x, y, 0.0, 1e-16, 1e-16, layout=nir.layout)
    _match(py, ft)


def test_invalid_inputs(twoband):
    _kernel, nir, art = twoband
    x, y = _binding(1, variant="x"), _binding(1, variant="y")
    cases = [
        (-0.1, 0.5, 0.0, 4),
        (1.1, 0.5, 0.0, 4),
        (0.5, 0.5, float("nan"), 3),
        (0.5, 0.5, float("inf"), 3),
    ]
    for u1, u2, log_q, code in cases:
        py = interpret_native_kernel(nir, x=x, y=y, log_q=log_q, u1=u1, u2=u2)
        ft = art.run_diagram(x, y, log_q, u1, u2, layout=nir.layout)
        assert py.status == code
        assert ft.status == code
        assert ft.accepted is False
    zero = Binding(k=x.k, q=x.q, tau=x.tau, t=x.t, omega=x.omega, g=0.0, L=x.L, delta=x.delta, gap=x.gap)
    py = interpret_native_kernel(nir, x=zero, y=y, log_q=0.0, u1=0.5, u2=0.5)
    ft = art.run_diagram(zero, y, 0.0, 0.5, 0.5, layout=nir.layout)
    assert py.status == 2
    assert ft.status == 2


def test_symmetric_proposal_rejects_nonzero_q(tmp_path):
    from fepdmc_toolkit.diagram_compiler.proposal import ProposalSpec

    kernel = compile_delayed_acceptance(build_scalar_ir(1), proposal_spec=ProposalSpec.symmetric())
    nir = lower_native_kernel(kernel)
    art = compile_native_kernel(nir, workdir=tmp_path / "sym")
    x, y = _binding(1, variant="x"), _binding(1, variant="y")
    py0 = interpret_native_kernel(nir, x=x, y=y, log_q=0.0, u1=1e-16, u2=1e-16)
    ft0 = art.run_diagram(x, y, 0.0, 1e-16, 1e-16, layout=nir.layout)
    _match(py0, ft0)
    py = interpret_native_kernel(nir, x=x, y=y, log_q=0.4, u1=1e-16, u2=1e-16)
    ft = art.run_diagram(x, y, 0.4, 1e-16, 1e-16, layout=nir.layout)
    assert py.status == 3 and ft.status == 3


def test_finite_state_from_compiled_score(score_art):
    nir, art = score_art
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
            ft = art.run_score(
                math.log(cheap[i]), math.log(cheap[j]), math.log(pi[i]), math.log(pi[j]), log_q, 1e-16, 1e-16
            )
            assert ft.status == 0
            py = interpret_native_kernel(
                nir,
                log_wx=math.log(cheap[i]),
                log_wy=math.log(cheap[j]),
                log_px=math.log(pi[i]),
                log_py=math.log(pi[j]),
                log_q=log_q,
                u1=1e-16,
                u2=1e-16,
            )
            _match(py, ft)
            a1 = 1.0 if ft.ell_hat >= 0.0 else math.exp(ft.ell_hat)
            delta = ft.ell_R - ft.ell_hat
            a2 = 1.0 if delta >= 0.0 else math.exp(delta)
            p[i, j] = q[i, j] * a1 * a2
        p[i, i] = 1.0 - p[i].sum()
    residual = max(abs(pi[i] * p[i, j] - pi[j] * p[j, i]) for i in range(n) for j in range(n))
    assert residual <= ATOL
    assert np.allclose(pi @ p, pi, atol=ATOL, rtol=RTOL)


def test_mutation_real32_kind_detected(twoband):
    _kernel, nir, _art = twoband
    source = emit_fortran(nir)
    assert "real(real64)" in source
    assert "real(real32)" not in source
    mutated = source.replace("real64", "real32")
    assert "real(real32)" in mutated
    x, y = _binding(1, variant="x"), _binding(1, variant="y")
    try:
        art = compile_fortran_source(mutated, driver=emit_driver(nir), workdir=None)
    except RuntimeError:
        return
    py = interpret_native_kernel(nir, x=x, y=y, log_q=0.0, u1=1e-16, u2=1e-16)
    ft = art.run_diagram(x, y, 0.0, 1e-16, 1e-16, layout=nir.layout)
    mismatch = (
        ft.accepted != py.accepted
        or abs(ft.ell_hat - py.ell_hat) > ATOL
        or (py.ell_R is not None and abs(ft.ell_R - py.ell_R) > ATOL)
    )
    assert mismatch


def test_mutation_stage2_drops_ell_hat(twoband, tmp_path):
    kernel, nir, _art = twoband
    source = emit_fortran(nir)
    mutated = source.replace("delta = ell_R - ell_hat", "delta = ell_R")
    assert mutated != source
    art = compile_fortran_source(mutated, driver=emit_driver(nir), workdir=tmp_path / "mut_delta")
    x, y = _binding(1, variant="x"), _binding(1, variant="y")
    py = interpret_native_kernel(nir, x=x, y=y, log_q=0.0, u1=1e-16, u2=0.8)
    ft = art.run_diagram(x, y, 0.0, 1e-16, 0.8, layout=nir.layout)
    assert py.accepted is True
    assert ft.accepted is False


def test_mutation_proposal_sign_and_field_order(twoband, tmp_path):
    _kernel, nir, _art = twoband
    source = emit_fortran(nir)
    mutated = source.replace("ell_R = diff_pi + log_q", "ell_R = diff_pi - log_q")
    assert mutated != source
    art = compile_fortran_source(mutated, driver=emit_driver(nir), workdir=tmp_path / "mut_q")
    x, y = _binding(1, variant="x"), _binding(1, variant="y")
    py = interpret_native_kernel(nir, x=x, y=y, log_q=0.4, u1=1e-16, u2=1e-16)
    ft = art.run_diagram(x, y, 0.4, 1e-16, 1e-16, layout=nir.layout)
    assert abs(ft.ell_R - py.ell_R) > ATOL
    kernel_mut = source.replace("x_k, x_q, x_tau, x_t", "x_t, x_q, x_tau, x_k", 1)
    assert kernel_mut != source
    art2 = compile_fortran_source(kernel_mut, driver=emit_driver(nir), workdir=tmp_path / "mut_abi")
    py2 = interpret_native_kernel(nir, x=x, y=y, log_q=0.0, u1=1e-16, u2=1e-16)
    ft2 = art2.run_diagram(x, y, 0.0, 1e-16, 1e-16, layout=nir.layout)
    assert abs(ft2.ell_hat - py2.ell_hat) > ATOL or ft2.accepted != py2.accepted or ft2.status != py2.status


def test_mutation_hoist_exact_and_status_and_target(twoband, tmp_path):
    from fepdmc_toolkit.diagram_compiler.fortran_compile import RUNTIME_FILE

    _kernel, nir, _art = twoband
    source = emit_fortran(nir)
    hoist = source.replace("call eval_graph_cheap(x_k,", "call eval_graph_exact(x_k,", 1)
    assert hoist != source
    with pytest.raises(NativeIRError):
        validate_generated_fortran(hoist, nir)
    art = compile_fortran_source(hoist, driver=emit_driver(nir), workdir=tmp_path / "mut_hoist")
    x, y = _binding(1, variant="x"), _binding(1, variant="y")
    ft = art.run_diagram(x, y, 0.0, 0.9, 0.1, layout=nir.layout)
    assert ft.reject_stage == 1
    assert ft.n_exact_graph > 0
    runtime = RUNTIME_FILE.read_text()
    status_rt = tmp_path / "runtime_status.f90"
    status_rt.write_text(runtime.replace("NK_INVALID_RANDOM_UNIFORM = 4_int64", "NK_INVALID_RANDOM_UNIFORM = 6_int64"))
    art_s = compile_fortran_source(source, driver=emit_driver(nir), runtime=status_rt, workdir=tmp_path / "mut_status")
    py = interpret_native_kernel(nir, x=x, y=y, log_q=0.0, u1=-0.1, u2=0.5)
    ft_s = art_s.run_diagram(x, y, 0.0, -0.1, 0.5, layout=nir.layout)
    assert py.status == 4
    assert ft_s.status == 6
    target_rt = tmp_path / "runtime_target.f90"
    target_rt.write_text(runtime.replace("re <= 0.0_real64", "re < -1.0e30_real64"))
    art_t = compile_fortran_source(source, driver=emit_driver(nir), runtime=target_rt, workdir=tmp_path / "mut_target")
    zero = Binding(k=x.k, q=x.q, tau=x.tau, t=x.t, omega=x.omega, g=0.0, L=x.L, delta=x.delta, gap=x.gap)
    py_z = interpret_native_kernel(nir, x=zero, y=y, log_q=0.0, u1=0.5, u2=0.5)
    ft_z = art_t.run_diagram(zero, y, 0.0, 0.5, 0.5, layout=nir.layout)
    assert py_z.status == 2
    assert ft_z.status != py_z.status


def test_o2_matches_o0(twoband, tmp_path):
    _kernel, nir, art0 = twoband
    art2 = compile_native_kernel(
        nir,
        workdir=tmp_path / "o2",
        flags=("-O2", "-std=f2008", "-ffree-line-length-none"),
    )
    x, y = _binding(1, variant="x"), _binding(1, variant="y")
    for log_q, u1, u2 in ((0.0, 0.9, 0.1), (0.0, 1e-16, 1e-16)):
        a = art0.run_diagram(x, y, log_q, u1, u2, layout=nir.layout)
        b = art2.run_diagram(x, y, log_q, u1, u2, layout=nir.layout)
        assert a.status == b.status
        assert a.accepted == b.accepted
        assert a.reject_stage == b.reject_stage
        assert a.exact_y_evaluated == b.exact_y_evaluated
        assert a.ell_hat == pytest.approx(b.ell_hat, abs=ATOL, rel=RTOL)


def test_native_ir_digest_recorded(twoband):
    _kernel, nir, _art = twoband
    source = emit_fortran(nir)
    assert native_ir_digest(nir) in source
    assert BACKEND in source
    assert "native_kernel_v1" in source
    assert "caller_guarantees_exact_x_valid_belongs_to_x" in source
