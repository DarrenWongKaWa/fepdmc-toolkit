"""Compile DAG Fortran plus the shared leaf-table evaluator. One .so per spec+tag."""

from __future__ import annotations

import ctypes
import hashlib
import shutil
import subprocess
import time
from pathlib import Path

import numpy as np

from fepdmc_toolkit.scheme_d.fortran_ir import emit_fortran

from .leaf_fortran import emit_leaf_assignments
from .leaves import LeafTable
from .plan_spec import PlanSpec


def fuse_module(
    dag_source: str,
    leaf_lines: list[str],
    batch_leaf_lines: list[str],
    spec: PlanSpec,
    n_input: int,
    n_work: int,
    n_output: int,
) -> str:
    body = dag_source.rstrip()
    if not body.endswith("end module"):
        raise ValueError("emit_fortran must end with end module")
    body = body[: -len("end module")].rstrip()
    n = spec.n
    ntau = 2 * n
    nw = max(1, n_work)
    nin = max(1, n_input)
    nout = max(1, n_output)
    phys = [
        'subroutine evaluate_physics(k, q, tau, t, omega, g, Lval, w, y, status) bind(C, name="scheme_d_eval_phys")',
        "real(c_double), intent(in), value :: k, t, omega, g, Lval",
        f"real(c_double), intent(in) :: q({n})",
        f"real(c_double), intent(in) :: tau({ntau})",
        f"complex(c_double_complex), intent(inout) :: w({nw})",
        f"complex(c_double_complex), intent(out) :: y({nout})",
        "integer(c_int), intent(out) :: status",
        f"complex(c_double_complex) :: x({nin})",
        "real(c_double) :: k_eff, dtau, xi, gv",
        "status=0",
        *leaf_lines,
        "call evaluate(x, w, y, status)",
        "end subroutine",
        'subroutine evaluate_physics_batch(nsamples, k_all, q_all, tau_all, t_all, omega_all, g_all, L_all, w, y, status) bind(C, name="scheme_d_eval_phys_batch")',
        "integer(c_int), intent(in), value :: nsamples",
        "real(c_double), intent(in) :: k_all(nsamples)",
        f"real(c_double), intent(in) :: q_all({n}, nsamples)",
        f"real(c_double), intent(in) :: tau_all({ntau}, nsamples)",
        "real(c_double), intent(in) :: t_all(nsamples), omega_all(nsamples), g_all(nsamples), L_all(nsamples)",
        f"complex(c_double_complex), intent(inout) :: w({nw})",
        "complex(c_double_complex), intent(out) :: y(nsamples)",
        "integer(c_int), intent(out) :: status(nsamples)",
        f"complex(c_double_complex) :: x({nin}), y1({nout})",
        f"real(c_double) :: k, t, omega, g, Lval, q({n}), tau({ntau})",
        "real(c_double) :: k_eff, dtau, xi, gv",
        "integer :: is, st",
        "do is=1, nsamples",
        "k = k_all(is)",
        "t = t_all(is)",
        "omega = omega_all(is)",
        "g = g_all(is)",
        "Lval = L_all(is)",
        f"q(1:{n}) = q_all(1:{n}, is)",
        f"tau(1:{ntau}) = tau_all(1:{ntau}, is)",
        "st = 0",
        *batch_leaf_lines,
        "call evaluate(x, w, y1, st)",
        "status(is) = st",
        "if (st == 0) then",
        "y(is) = y1(1)",
        "else",
        "y(is) = cmplx(ieee_value(0.0d0, ieee_quiet_nan), 0.0d0, kind=c_double)",
        "end if",
        "end do",
        "end subroutine",
        "end module",
    ]
    return body + "\n" + "\n".join(phys) + "\n"


class FusedContext:
    def __init__(self, kernel: FusedKernel) -> None:
        self.kernel = kernel
        self.kernel_id = id(kernel)
        self.spec_id = kernel.spec.spec_id
        self.constructor = kernel.tag
        self.work = np.empty(max(1, kernel.ir.n_work), dtype=np.complex128)
        self.output = np.empty(max(1, kernel.ir.n_output), dtype=np.complex128)
        self.status = ctypes.c_int()
        self.q = np.empty(kernel.spec.n, dtype=np.float64)
        self.tau = np.empty(2 * kernel.spec.n, dtype=np.float64)

    def fill_sample(self, binding) -> None:
        self.q[:] = binding.q
        self.tau[:] = binding.tau

    def run_filled(self, binding) -> np.ndarray:
        self.kernel.fn_phys(
            ctypes.c_double(binding.k),
            self.q.ctypes.data,
            self.tau.ctypes.data,
            ctypes.c_double(binding.t),
            ctypes.c_double(binding.omega),
            ctypes.c_double(binding.g),
            ctypes.c_double(float(binding.L)),
            self.work.ctypes.data,
            self.output.ctypes.data,
            ctypes.byref(self.status),
        )
        if self.status.value:
            raise FloatingPointError(f"fused evaluator status={self.status.value}")
        return self.output.copy()

    def run(self, binding) -> np.ndarray:
        self.fill_sample(binding)
        return self.run_filled(binding)


class FusedKernel:
    def __init__(
        self, spec: PlanSpec, graph, leaves: LeafTable, directory: Path, compiler=None, *, tag: str = "R1-symbolic"
    ):
        self.spec = spec
        self.tag = tag
        self.graph = graph
        self.leaves = leaves
        self.ir = emit_fortran(graph)
        leaf_lines = emit_leaf_assignments(leaves, spec, self.ir.variables, include_l_guard=True)
        batch_leaf_lines = emit_leaf_assignments(leaves, spec, self.ir.variables, include_l_guard=False)
        self.leaf_fortran = "\n".join(leaf_lines) + "\n"
        self.source = fuse_module(
            self.ir.source,
            leaf_lines,
            batch_leaf_lines,
            spec,
            self.ir.n_input,
            self.ir.n_work,
            self.ir.n_output,
        )
        compiler = compiler or shutil.which("gfortran")
        if not compiler:
            raise FileNotFoundError("gfortran unavailable")
        directory = Path(directory)
        directory.mkdir(parents=True, exist_ok=True)
        digest = hashlib.sha256((self.source + str(compiler) + " -O2 -fno-fast-math fused").encode()).hexdigest()[:20]
        self.source_path = directory / f"fused_{digest}.f90"
        self.library_path = directory / f"fused_{digest}.so"
        self.source_path.write_text(self.source)
        start = time.perf_counter()
        proc = subprocess.run(
            [
                compiler,
                "-shared",
                "-fPIC",
                "-O2",
                "-fno-fast-math",
                "-ffree-line-length-none",
                "-J",
                str(directory),
                str(self.source_path),
                "-o",
                str(self.library_path),
            ],
            text=True,
            capture_output=True,
            timeout=180,
        )
        self.compile_s = time.perf_counter() - start
        if proc.returncode:
            raise RuntimeError(proc.stderr)
        self.lib = ctypes.CDLL(str(self.library_path.resolve()))
        self.fn_phys = self.lib.scheme_d_eval_phys
        self.fn_phys.argtypes = [
            ctypes.c_double,
            ctypes.c_void_p,
            ctypes.c_void_p,
            ctypes.c_double,
            ctypes.c_double,
            ctypes.c_double,
            ctypes.c_double,
            ctypes.c_void_p,
            ctypes.c_void_p,
            ctypes.POINTER(ctypes.c_int),
        ]
        self.fn_phys.restype = None
        self.fn_batch = self.lib.scheme_d_eval_phys_batch
        self.fn_batch.argtypes = [
            ctypes.c_int,
            ctypes.c_void_p,
            ctypes.c_void_p,
            ctypes.c_void_p,
            ctypes.c_void_p,
            ctypes.c_void_p,
            ctypes.c_void_p,
            ctypes.c_void_p,
            ctypes.c_void_p,
            ctypes.c_void_p,
            ctypes.c_void_p,
        ]
        self.fn_batch.restype = None

    def spawn(self) -> FusedContext:
        return FusedContext(self)

    def run_batch(self, k, q, tau, t, omega, g, L, work, y, status) -> None:
        n = int(k.shape[0])
        if n < 1:
            raise ValueError("nsamples>=1")
        self.fn_batch(
            ctypes.c_int(n),
            k.ctypes.data,
            q.ctypes.data,
            tau.ctypes.data,
            t.ctypes.data,
            omega.ctypes.data,
            g.ctypes.data,
            L.ctypes.data,
            work.ctypes.data,
            y.ctypes.data,
            status.ctypes.data,
        )
