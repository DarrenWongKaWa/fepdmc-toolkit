"""Stage 5: the Fortran backend, compiled with the local gfortran.

The two-band n = 1 kernel is emitted as Fortran, compiled together with the
primitive runtime and run on every path. The discrete results must equal those of
the NativeKernelIR interpreter and the scores must agree within the tolerance; each
runtime primitive must agree with its Python reference; the compiled score kernel
must keep detailed balance on the three-state chain; the emitted source must be
deterministic and equal the validated one (validation/diagram_compiler/golden.json).

Records: generated_twoband_n1.f90, kernel_differential.csv, call_counts.csv,
primitive_differential.csv, finite_state_fortran.csv, source_provenance.json,
FORTRAN_BACKEND_SPEC.json and summary.json. Needs gfortran (FEPDMC_TOOLKIT_GFORTRAN
selects another compiler).
"""

from __future__ import annotations

import argparse
import hashlib
import math
import sys
import tempfile
from pathlib import Path

import numpy as np
from _common import (
    BALANCE_TOL,
    CHEAP,
    PI,
    Q_ASYMMETRIC,
    ROOT,
    balance_residuals,
    finish,
    golden,
    twoband_binding,
    write_csv,
    write_json,
)

from fepdmc_toolkit.diagram_compiler import (
    build_twoband_ir,
    compile_delayed_acceptance,
    emit_fortran,
    interpret_native_kernel,
    lower_native_kernel,
    lower_native_score_kernel,
)
from fepdmc_toolkit.diagram_compiler.fortran_codegen import BACKEND, native_ir_digest
from fepdmc_toolkit.diagram_compiler.fortran_compile import (
    ATOL,
    RTOL,
    RUNTIME_FILE,
    compile_native_kernel,
    compile_primitive_driver,
    compiler_version,
    find_gfortran,
    gfortran_available,
    run_primitive,
)
from fepdmc_toolkit.diagram_compiler.native_interpret import _eval_op
from fepdmc_toolkit.diagram_compiler.native_ir import TYPE_C128, TYPE_C128_M2, NativeOp


def _parse_z(stdout):
    re, im = map(float, stdout.splitlines()[0].split())
    return complex(re, im)


def _parse_m(stdout):
    vals = list(map(float, stdout.splitlines()[0].split()))
    return np.array(
        [[vals[0] + 1j * vals[1], vals[2] + 1j * vals[3]], [vals[4] + 1j * vals[5], vals[6] + 1j * vals[7]]],
        dtype=np.complex128,
    )


def _dump_m(matrix):
    return " ".join(f"{matrix[i, j].real} {matrix[i, j].imag}" for i in range(2) for j in range(2))


def _within(err, ref) -> bool:
    return err is None or err <= ATOL + RTOL * abs(ref)


def main(output: Path) -> int:
    output.mkdir(parents=True, exist_ok=True)
    kernel = compile_delayed_acceptance(build_twoband_ir(1))
    nir = lower_native_kernel(kernel)
    source = emit_fortran(nir)
    (output / "generated_twoband_n1.f90").write_text(source, encoding="utf-8", newline="\n")
    art = compile_native_kernel(nir, workdir=tempfile.mkdtemp(prefix="stage5_tb_"))
    x, y = twoband_binding("x"), twoband_binding("y")
    probe = kernel.evaluate_transition(x, y, 0.0, 1e-16, 1e-16)
    alpha2 = math.exp(min(0.0, probe.ell_R - probe.ell_hat))
    cases = [
        ("stage1_reject", 0.0, 0.9, 0.1),
        ("stage2_reject", 0.0, 1e-16, (alpha2 + 1.0) / 2.0),
        ("accept", 0.0, 1e-16, 1e-16),
        ("asymmetric_q", 0.4, 1e-16, 1e-16),
    ]
    kernel_rows = []
    call_rows = []
    scores_ok = True
    for name, log_q, u1, u2 in cases:
        py = interpret_native_kernel(nir, x=x, y=y, log_q=log_q, u1=u1, u2=u2)
        ft = art.run_diagram(x, y, log_q, u1, u2, layout=nir.layout)
        hat_err = abs(ft.ell_hat - py.ell_hat) if py.status == 0 else None
        r_err = abs(ft.ell_R - py.ell_R) if py.ell_R is not None else None
        scores_ok = scores_ok and _within(hat_err, py.ell_hat if py.status == 0 else 0.0)
        scores_ok = scores_ok and _within(r_err, py.ell_R if py.ell_R is not None else 0.0)
        kernel_rows.append(
            {
                "name": name,
                "status_match": py.status == ft.status,
                "accepted_match": py.accepted == ft.accepted,
                "reject_stage_match": py.reject_stage == ft.reject_stage,
                "exact_y_match": py.exact_y_evaluated == ft.exact_y_evaluated,
                "ell_hat_abs_err": hat_err,
                "ell_R_abs_err": r_err,
                "py_accepted": py.accepted,
                "ft_accepted": ft.accepted,
                "py_reject_stage": py.reject_stage,
                "ft_reject_stage": ft.reject_stage,
                "py_ell_hat": py.ell_hat,
                "ft_ell_hat": ft.ell_hat,
                "py_ell_R": py.ell_R,
                "ft_ell_R": ft.ell_R if ft.ell_R_valid else None,
                "py_exact_y": py.exact_y_evaluated,
                "ft_exact_y": ft.exact_y_evaluated,
            }
        )
        call_rows.append(
            {
                "name": name,
                "n_exact_graph": ft.n_exact_graph,
                "n_electron_twoband": ft.n_electron_twoband,
                "n_vertex_sigmaz": ft.n_vertex_sigmaz,
                "n_matmul": ft.n_matmul,
                "n_trace": ft.n_trace,
                "stage1_reject": name == "stage1_reject",
            }
        )
    write_csv(output / "kernel_differential.csv", kernel_rows)
    write_csv(output / "call_counts.csv", call_rows)

    prim = compile_primitive_driver(workdir=tempfile.mkdtemp(prefix="stage5_prim_"))
    prim_rows = []

    def add_z(name, stdin, ref):
        got = _parse_z(run_primitive(prim, stdin))
        prim_rows.append(
            {
                "name": name,
                "abs_err": abs(got - ref),
                "within_tol": abs(got - ref) <= ATOL + RTOL * abs(ref),
            }
        )

    ref = _eval_op(
        NativeOp("z", "prim.electron_scalar", TYPE_C128, ("a", "b", "c")),
        {"a": 0.3, "b": 1.0, "c": 0.05},
        binding=None,
        inputs={},
    )
    add_z("electron_scalar", "1\n0.3 1.0 0.05\n", ref)
    ref = _eval_op(
        NativeOp("z", "prim.phonon", TYPE_C128, ("a", "b")),
        {"a": 0.8, "b": 0.05},
        binding=None,
        inputs={},
    )
    add_z("phonon", "3\n0.8 0.05\n", ref)
    ref = _eval_op(
        NativeOp("z", "prim.prefactor", TYPE_C128, ("g", "L"), {"n": 1}),
        {"g": 0.5, "L": 4.0},
        binding=None,
        inputs={},
    )
    add_z("prefactor", "4\n0.5 4.0 1\n", ref)
    u_ref = _eval_op(
        NativeOp("z", "prim.electron_twoband", TYPE_C128_M2, ("a", "b", "c", "d", "e")),
        {"a": 0.3, "b": 1.0, "c": 0.05, "d": 0.35, "e": 0.5},
        binding=None,
        inputs={},
    )
    u_got = _parse_m(run_primitive(prim, "2\n0.3 1.0 0.05 0.35 0.5\n"))
    prim_rows.append(
        {
            "name": "electron_twoband",
            "abs_err": float(np.max(np.abs(u_got - u_ref))),
            "within_tol": bool(np.allclose(u_got, u_ref, atol=ATOL, rtol=RTOL)),
        }
    )
    v_ref = _eval_op(
        NativeOp("z", "prim.vertex_sigmaz", TYPE_C128_M2, ("g", "L")),
        {"g": 0.5, "L": 4.0},
        binding=None,
        inputs={},
    )
    v_got = _parse_m(run_primitive(prim, "5\n0.5 4.0\n"))
    prim_rows.append(
        {
            "name": "vertex_sigmaz",
            "abs_err": float(np.max(np.abs(v_got - v_ref))),
            "within_tol": bool(np.allclose(v_got, v_ref, atol=ATOL, rtol=RTOL)),
        }
    )
    a = np.array([[1 + 2j, 3 - 4j], [5 + 0j, -1 + 7j]], dtype=np.complex128)
    b = np.array([[0.2 - 0.1j, 1j], [2 + 3j, 4 - 5j]], dtype=np.complex128)
    m_got = _parse_m(run_primitive(prim, f"6\n{_dump_m(a)}\n{_dump_m(b)}\n"))
    prim_rows.append(
        {
            "name": "matmul",
            "abs_err": float(np.max(np.abs(m_got - (a @ b)))),
            "within_tol": bool(np.allclose(m_got, a @ b, atol=ATOL, rtol=RTOL)),
        }
    )
    s = 0.3 - 0.4j
    s_got = _parse_m(run_primitive(prim, f"7\n{_dump_m(a)}\n{s.real} {s.imag}\n"))
    prim_rows.append(
        {
            "name": "scale_m2",
            "abs_err": float(np.max(np.abs(s_got - a * s))),
            "within_tol": bool(np.allclose(s_got, a * s, atol=ATOL, rtol=RTOL)),
        }
    )
    tr = _parse_z(run_primitive(prim, f"8\n{_dump_m(a)}\n"))
    prim_rows.append(
        {
            "name": "trace",
            "abs_err": abs(tr - complex(np.trace(a))),
            "within_tol": abs(tr - complex(np.trace(a))) <= ATOL,
        }
    )
    write_csv(output / "primitive_differential.csv", prim_rows)

    score = lower_native_score_kernel()
    score_art = compile_native_kernel(score, workdir=tempfile.mkdtemp(prefix="stage5_score_"))
    n = len(PI)
    p = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            log_q = math.log(Q_ASYMMETRIC[j, i]) - math.log(Q_ASYMMETRIC[i, j])
            ft = score_art.run_score(
                math.log(CHEAP[i]), math.log(CHEAP[j]), math.log(PI[i]), math.log(PI[j]), log_q, 1e-16, 1e-16
            )
            a1 = 1.0 if ft.ell_hat >= 0.0 else math.exp(ft.ell_hat)
            delta = ft.ell_R - ft.ell_hat
            a2 = 1.0 if delta >= 0.0 else math.exp(delta)
            p[i, j] = Q_ASYMMETRIC[i, j] * a1 * a2
        p[i, i] = 1.0 - p[i].sum()
    residual, stat = balance_residuals(p)
    write_csv(
        output / "finite_state_fortran.csv",
        [{"balance_residual": residual, "stationarity_residual": stat, "atol": ATOL, "rtol": RTOL}],
    )

    compiler = find_gfortran()
    flags = ("-O0", "-std=f2008", "-ffree-line-length-none", "-fcheck=bounds")
    source_sha = hashlib.sha256(source.encode("utf-8")).hexdigest()
    repeat_sha = hashlib.sha256(emit_fortran(nir).encode("utf-8")).hexdigest()
    runtime_sha = hashlib.sha256(RUNTIME_FILE.read_bytes()).hexdigest()
    write_json(
        output / "source_provenance.json",
        {
            "backend": BACKEND,
            "binary_sha256": art.binary_sha256,
            "compiler": compiler,
            "compiler_command": [
                compiler,
                *flags,
                "-J",
                ".",
                "native_kernel_runtime_v1.f90",
                "generated_da_kernel.f90",
                "nk_eval_driver.f90",
                "-o",
                "nk_eval",
            ],
            "compiler_version": compiler_version(compiler),
            "fortran_source_sha256": source_sha,
            "native_ir_digest": native_ir_digest(nir),
            "primitive_runtime": "native_kernel_v1",
            "runtime_sha256": runtime_sha,
            "repeat_source_sha256": repeat_sha,
        },
    )
    write_json(
        output / "FORTRAN_BACKEND_SPEC.json",
        {
            "abi_bool": "logical internally and on kernel dummy list; driver I/O uses integer 0/1",
            "structure": "generated kernel plus versioned primitive runtime",
            "backend": BACKEND,
            "cache_contract": "the caller guarantees that exact_x_valid belongs to x",
            "electron_twoband": "closed form for the 2x2 H = m I + K with K^2 = r^2 I (cosh, sinh); no LAPACK",
            "matrix_layout": "mathematical (i,j) stored as Fortran a(i,j) column-major",
            "primitive_set": "native_kernel_v1",
            "schema_version": 1,
            "status_codes": {
                "CACHE_KEY_MISMATCH": 5,
                "INVALID_CHEAP_WEIGHT": 2,
                "INVALID_EXACT_TARGET": 1,
                "INVALID_IR": 6,
                "INVALID_PROPOSAL_RATIO": 3,
                "INVALID_RANDOM_UNIFORM": 4,
                "OK": 0,
            },
            "tolerances": {"atol": ATOL, "discrete": "exact", "rtol": RTOL},
            "types": {
                "bool": "logical",
                "c128": "complex(real64)",
                "c128_m2": "complex(real64) :: a(2,2)",
                "f64": "real(real64)",
                "i64": "integer(int64)",
            },
        },
    )
    reference = golden()
    stage1 = next(r for r in call_rows if r["stage1_reject"])
    return finish(
        "stage 5 (Fortran backend)",
        output,
        {
            "compiled_kernel_matches_the_interpreter": all(
                r["status_match"] and r["accepted_match"] and r["reject_stage_match"] and r["exact_y_match"]
                for r in kernel_rows
            )
            and scores_ok,
            "no_exact_primitive_on_the_stage1_reject_path": stage1["n_exact_graph"] == 0
            and stage1["n_electron_twoband"] == 0
            and stage1["n_trace"] == 0,
            "every_primitive_matches_its_reference": all(r["within_tol"] for r in prim_rows),
            "detailed_balance_on_the_three_state_chain": residual <= BALANCE_TOL and stat <= BALANCE_TOL,
            "emitted_source_is_deterministic": source_sha == repeat_sha,
            "emitted_source_equals_the_validated_source": source_sha == reference["fortran_source_sha256"],
            "runtime_equals_the_validated_runtime": runtime_sha == reference["runtime_sha256"],
        },
        compiler_version=compiler_version(compiler),
        max_balance_residual=residual,
        max_primitive_abs_err=max(float(r["abs_err"]) for r in prim_rows),
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--output", type=Path, default=ROOT / "build/validation/stage5_fortran_backend")
    args = parser.parse_args()
    if not gfortran_available():
        parser.error("gfortran not found; install it or set FEPDMC_TOOLKIT_GFORTRAN")
    sys.exit(main(args.output))
