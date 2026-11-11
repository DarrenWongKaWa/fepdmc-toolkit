"""Stage 5: native_kernel_v1 Fortran primitives vs Python interpreter oracles."""

from __future__ import annotations

import math

import numpy as np
import pytest

from fepdmc_toolkit.diagram_compiler.fortran_compile import (
    ATOL,
    RTOL,
    compile_primitive_driver,
    gfortran_available,
    run_primitive,
)
from fepdmc_toolkit.diagram_compiler.native_interpret import _eval_op, _xi
from fepdmc_toolkit.diagram_compiler.native_ir import TYPE_C128, TYPE_C128_M2, NativeOp

pytestmark = pytest.mark.skipif(not gfortran_available(), reason="gfortran not found")


@pytest.fixture(scope="module")
def prim_exe(tmp_path_factory):
    return compile_primitive_driver(workdir=tmp_path_factory.mktemp("prim"))


def _c128(op_name, args, env):
    op = NativeOp("z", op_name, TYPE_C128, tuple(env))
    mapping = dict(zip(env, args))
    return _eval_op(op, mapping, binding=None, inputs={})


def _m2(op_name, arg_names, env):
    op = NativeOp("z", op_name, TYPE_C128_M2, tuple(arg_names))
    return _eval_op(op, env, binding=None, inputs={})


def _parse_z(stdout):
    re, im = map(float, stdout.splitlines()[0].split())
    return complex(re, im)


def _parse_m(stdout):
    vals = list(map(float, stdout.splitlines()[0].split()))
    assert len(vals) == 8
    return np.array(
        [[vals[0] + 1j * vals[1], vals[2] + 1j * vals[3]], [vals[4] + 1j * vals[5], vals[6] + 1j * vals[7]]],
        dtype=np.complex128,
    )


def _dump_m(matrix):
    return " ".join(f"{matrix[i, j].real} {matrix[i, j].imag}" for i in range(2) for j in range(2))


def _close_z(got, ref):
    assert got.real == pytest.approx(ref.real, abs=ATOL, rel=RTOL)
    assert got.imag == pytest.approx(ref.imag, abs=ATOL, rel=RTOL)


def _close_m(got, ref):
    assert np.allclose(got, ref, atol=ATOL, rtol=RTOL)


def test_electron_scalar_and_xi(prim_exe):
    cases = [(0.0, 1.0, 0.05), (0.3, 1.0, 0.05), (math.pi, 1.1, -0.07), (1.2, 0.5, 1.0)]
    for k, t, dtau in cases:
        ref = _c128("prim.electron_scalar", (k, t, dtau), ("a", "b", "c"))
        got = _parse_z(run_primitive(prim_exe, f"1\n{k} {t} {dtau}\n"))
        _close_z(got, ref)
        xi = float(run_primitive(prim_exe, f"13\n{k} {t}\n").splitlines()[0])
        assert xi == pytest.approx(_xi(k, t), abs=ATOL, rel=RTOL)


def test_electron_twoband_closed_form(prim_exe):
    rng = np.random.default_rng(5)
    grid = [
        (0.0, 1.0, 0.05, 0.35, 0.5),
        (0.3, 1.0, 0.05, 0.0, 0.5),
        (math.pi, 1.1, -0.07, 0.35, -0.8),
        (0.0, 1.0, 0.2, 0.0, 0.0),
    ]
    randoms = [
        (float(k), float(t), float(dt), float(d), float(g))
        for k, t, dt, d, g in rng.normal(size=(40, 5)) * np.array([1.5, 0.4, 0.2, 0.5, 0.7])
        + np.array([0.5, 1.0, 0.05, 0.2, 0.3])
    ]
    max_err = 0.0
    for k, t, dtau, delta, gap in grid + randoms:
        ref = _m2(
            "prim.electron_twoband",
            ("a", "b", "c", "d", "e"),
            {"a": k, "b": t, "c": dtau, "d": delta, "e": gap},
        )
        got = _parse_m(run_primitive(prim_exe, f"2\n{k} {t} {dtau} {delta} {gap}\n"))
        max_err = max(max_err, float(np.max(np.abs(got - ref))))
        _close_m(got, ref)
    assert max_err <= ATOL


def test_phonon_prefactor_vertex(prim_exe):
    ref = _c128("prim.phonon", (0.8, 0.05), ("a", "b"))
    got = _parse_z(run_primitive(prim_exe, "3\n0.8 0.05\n"))
    _close_z(got, ref)
    op = NativeOp("z", "prim.prefactor", TYPE_C128, ("g", "L"), {"n": 1})
    ref = _eval_op(op, {"g": 0.5, "L": 4.0}, binding=None, inputs={})
    got = _parse_z(run_primitive(prim_exe, "4\n0.5 4.0 1\n"))
    _close_z(got, ref)
    ref = _m2("prim.vertex_sigmaz", ("g", "L"), {"g": 0.5, "L": 4.0})
    got = _parse_m(run_primitive(prim_exe, "5\n0.5 4.0\n"))
    _close_m(got, ref)
    assert got[0, 0].real > 0 and got[1, 1].real < 0


def test_matmul_scale_trace_nonsymmetric(prim_exe):
    a = np.array([[1 + 2j, 3 - 4j], [5 + 0j, -1 + 7j]], dtype=np.complex128)
    b = np.array([[0.2 - 0.1j, 1j], [2 + 3j, 4 - 5j]], dtype=np.complex128)
    got = _parse_m(run_primitive(prim_exe, f"6\n{_dump_m(a)}\n{_dump_m(b)}\n"))
    _close_m(got, a @ b)
    assert not np.allclose(got, b @ a, atol=1e-8)
    assert not np.allclose(got, a.T @ b, atol=1e-8)
    s = 0.3 - 0.4j
    got_s = _parse_m(run_primitive(prim_exe, f"7\n{_dump_m(a)}\n{s.real} {s.imag}\n"))
    _close_m(got_s, a * s)
    tr = _parse_z(run_primitive(prim_exe, f"8\n{_dump_m(a)}\n"))
    _close_z(tr, complex(np.trace(a)))
    assert abs(tr - (a[0, 1] + a[1, 0])) > 1e-8
    got_add = _parse_m(run_primitive(prim_exe, f"9\n{_dump_m(a)}\n{_dump_m(b)}\n"))
    _close_m(got_add, a + b)


def test_log_positive_real_and_uniforms(prim_exe):
    out = run_primitive(prim_exe, "10\n2.0 0.0 1e-12\nINVALID_EXACT_TARGET\n")
    status, logw = out.splitlines()[0].split()
    assert int(status) == 0
    assert float(logw) == pytest.approx(math.log(2.0), abs=ATOL, rel=RTOL)
    out = run_primitive(prim_exe, "10\n-1.0 0.0 1e-12\nINVALID_CHEAP_WEIGHT\n")
    status, _logw = out.splitlines()[0].split()
    assert int(status) == 2
    out = run_primitive(prim_exe, "10\n1.0 1e-6 1e-12\nINVALID_EXACT_TARGET\n")
    status, _logw = out.splitlines()[0].split()
    assert int(status) == 1
    out = run_primitive(prim_exe, "11\n0.5 1e-300\n")
    status, logu = out.splitlines()[0].split()
    assert int(status) == 0
    assert float(logu) == pytest.approx(math.log(0.5), abs=ATOL, rel=RTOL)
    out = run_primitive(prim_exe, "11\n-0.1 1e-300\n")
    assert int(out.splitlines()[0].split()[0]) == 4
    out = run_primitive(prim_exe, "12\n0.4 1\n")
    assert int(out.splitlines()[0].split()[0]) == 3
    out = run_primitive(prim_exe, "12\n0.0 1\n")
    assert int(out.splitlines()[0].split()[0]) == 0


def test_mutation_runtime_matmul_transpose_and_trace(tmp_path):
    from fepdmc_toolkit.diagram_compiler.fortran_compile import RUNTIME_FILE, compile_primitive_driver, run_primitive

    runtime = RUNTIME_FILE.read_text()
    mutated = runtime.replace("c = matmul(a, b)", "c = matmul(b, a)")
    mutated = mutated.replace("z = a(1, 1) + a(2, 2)", "z = a(1, 2) + a(2, 1)")
    path = tmp_path / "runtime_mut.f90"
    path.write_text(mutated)
    exe = compile_primitive_driver(workdir=tmp_path / "prim_mut", runtime=path)
    a = np.array([[1 + 2j, 3 - 4j], [5 + 0j, -1 + 7j]], dtype=np.complex128)
    b = np.array([[0.2 - 0.1j, 1j], [2 + 3j, 4 - 5j]], dtype=np.complex128)
    got = _parse_m(run_primitive(exe, f"6\n{_dump_m(a)}\n{_dump_m(b)}\n"))
    assert np.allclose(got, b @ a, atol=ATOL, rtol=RTOL)
    assert not np.allclose(got, a @ b, atol=1e-8)
    tr = _parse_z(run_primitive(exe, f"8\n{_dump_m(a)}\n"))
    assert abs(tr - complex(np.trace(a))) > 1e-8
    assert tr == pytest.approx(complex(a[0, 1] + a[1, 0]), abs=ATOL, rel=RTOL)


def test_mutation_wrong_tau_sign_and_vertex_and_prefactor(prim_exe):
    k, t, dtau = 0.3, 1.0, 0.2
    ref = _parse_z(run_primitive(prim_exe, f"1\n{k} {t} {dtau}\n"))
    flipped = _parse_z(run_primitive(prim_exe, f"1\n{k} {t} {-dtau}\n"))
    assert abs(ref - flipped) > 1e-8
    v = _parse_m(run_primitive(prim_exe, "5\n0.5 4.0\n"))
    assert abs(v[0, 0] + v[1, 1]) <= ATOL
    n1 = _parse_z(run_primitive(prim_exe, "4\n0.5 4.0 1\n"))
    n2 = _parse_z(run_primitive(prim_exe, "4\n0.5 4.0 2\n"))
    assert abs(n2 - n1 * n1) <= 1e-12 or abs(n2 - (n1 * n1)) / abs(n1 * n1) <= 1e-10
