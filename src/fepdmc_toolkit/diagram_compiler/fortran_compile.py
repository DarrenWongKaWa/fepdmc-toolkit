"""Compile generated NativeKernelIR Fortran. Test/script helper, not FEP-DMC."""

from __future__ import annotations

import hashlib
import os
import shutil
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path

from .evaluator import Binding
from .fortran_codegen import BACKEND, emit_driver, emit_fortran, native_ir_digest
from .native_ir import NativeKernelIR

RUNTIME_DIR = Path(__file__).resolve().parent / "fortran"
RUNTIME_FILE = RUNTIME_DIR / "native_kernel_runtime_v1.f90"
PRIM_DRIVER = RUNTIME_DIR / "nk_prim_driver.f90"

DEFAULT_FLAGS = ("-O0", "-std=f2008", "-ffree-line-length-none", "-fcheck=bounds")
ATOL = 1e-12
RTOL = 1e-10


def find_gfortran() -> str:
    explicit = os.environ.get("FEPDMC_TOOLKIT_GFORTRAN")
    if explicit:
        return explicit
    path = shutil.which("gfortran")
    if path:
        return path
    homebrew = Path("/opt/homebrew/bin/gfortran")
    if homebrew.exists():
        return str(homebrew)
    raise FileNotFoundError("gfortran not found")


def gfortran_available() -> bool:
    """True if a Fortran compiler is found (FEPDMC_TOOLKIT_GFORTRAN, gfortran on PATH or Homebrew)."""
    try:
        find_gfortran()
    except FileNotFoundError:
        return False
    return True


def compiler_version(compiler: str | None = None) -> str:
    path = compiler or find_gfortran()
    result = subprocess.run([path, "--version"], check=True, capture_output=True, text=True)
    return result.stdout.splitlines()[0]


@dataclass(frozen=True)
class FortranResult:
    status: int
    accepted: bool
    reject_stage: int
    ell_hat: float
    ell_R_valid: bool
    ell_R: float
    exact_y_evaluated: bool
    n_electron_twoband: int
    n_vertex_sigmaz: int
    n_matmul: int
    n_trace: int
    n_exact_graph: int
    stdout: str


@dataclass
class CompiledFortran:
    exe: Path
    source_path: Path
    driver_path: Path
    workdir: Path
    compiler: str
    flags: tuple[str, ...]
    source_sha256: str
    binary_sha256: str
    nir_digest: str
    kind: str

    def run_text(self, stdin: str) -> FortranResult:
        proc = subprocess.run(
            [str(self.exe)],
            input=stdin,
            text=True,
            capture_output=True,
            check=False,
            cwd=str(self.workdir),
        )
        if proc.returncode != 0:
            raise RuntimeError(f"native kernel exited {proc.returncode}: {proc.stderr}\nstdin={stdin!r}")
        line = proc.stdout.strip().splitlines()[-1]
        parts = line.split()
        if len(parts) != 12:
            raise RuntimeError(f"unexpected driver output {proc.stdout!r}")
        status = int(parts[0])
        return FortranResult(
            status=status,
            accepted=bool(int(parts[1])),
            reject_stage=int(parts[2]),
            ell_hat=float(parts[3]),
            ell_R_valid=bool(int(parts[4])),
            ell_R=float(parts[5]),
            exact_y_evaluated=bool(int(parts[6])),
            n_electron_twoband=int(parts[7]),
            n_vertex_sigmaz=int(parts[8]),
            n_matmul=int(parts[9]),
            n_trace=int(parts[10]),
            n_exact_graph=int(parts[11]),
            stdout=proc.stdout,
        )

    def run_diagram(
        self,
        x: Binding,
        y: Binding,
        log_q: float,
        u1: float,
        u2: float,
        *,
        exact_x_valid: bool = False,
        exact_log_weight_x: float = 0.0,
        layout=None,
    ) -> FortranResult:
        if layout is None:
            raise ValueError("layout required")
        chunks: list[str] = []
        for binding in (x, y):
            for field in layout.fields:
                value = getattr(binding, field.name)
                if field.shape == ():
                    chunks.append(str(value))
                else:
                    chunks.append(" ".join(str(float(v)) for v in value))
        chunks.append(f"{float(log_q)} {float(u1)} {float(u2)}")
        chunks.append(f"{1 if exact_x_valid else 0} {float(exact_log_weight_x)}")
        return self.run_text("\n".join(chunks) + "\n")

    def run_score(self, log_wx, log_wy, log_px, log_py, log_q, u1, u2) -> FortranResult:
        stdin = (
            f"{float(log_wx)} {float(log_wy)} {float(log_px)} {float(log_py)}\n{float(log_q)} {float(u1)} {float(u2)}\n"
        )
        return self.run_text(stdin)


def compile_fortran_source(
    source: str,
    *,
    runtime: Path | str | None = None,
    driver: str | None = None,
    workdir: Path | str | None = None,
    compiler: str | None = None,
    flags: tuple[str, ...] = DEFAULT_FLAGS,
    exe_name: str = "nk_eval",
) -> CompiledFortran:
    compiler = compiler or find_gfortran()
    runtime_path = Path(runtime) if runtime is not None else RUNTIME_FILE
    if driver is None:
        raise ValueError("driver source is required")
    if workdir is None:
        workdir = Path(tempfile.mkdtemp(prefix="fepdmc_toolkit_"))
    else:
        workdir = Path(workdir)
        workdir.mkdir(parents=True, exist_ok=True)
    src_path = workdir / "generated_da_kernel.f90"
    drv_path = workdir / "nk_eval_driver.f90"
    rt_copy = workdir / "native_kernel_runtime_v1.f90"
    src_path.write_text(source, encoding="utf-8", newline="\n")
    drv_path.write_text(driver, encoding="utf-8", newline="\n")
    shutil.copyfile(runtime_path, rt_copy)
    exe = workdir / exe_name
    cmd = [compiler, *flags, "-J", str(workdir), str(rt_copy), str(src_path), str(drv_path), "-o", str(exe)]
    proc = subprocess.run(cmd, capture_output=True, text=True, check=False)
    if proc.returncode != 0:
        raise RuntimeError(f"gfortran failed: {proc.stderr}\ncommand: {cmd}")
    source_sha = hashlib.sha256(source.encode("utf-8")).hexdigest()
    binary_sha = hashlib.sha256(exe.read_bytes()).hexdigest()
    digest = ""
    for line in source.splitlines():
        if "native_ir_digest=" in line:
            digest = line.split("native_ir_digest=", 1)[1].strip()
            break
    kind = "score" if "log_wx_in" in source else "diagram"
    return CompiledFortran(
        exe=exe,
        source_path=src_path,
        driver_path=drv_path,
        workdir=workdir,
        compiler=compiler,
        flags=tuple(flags),
        source_sha256=source_sha,
        binary_sha256=binary_sha,
        nir_digest=digest,
        kind=kind,
    )


def compile_native_kernel(
    nir: NativeKernelIR,
    *,
    backend: str = BACKEND,
    workdir: Path | str | None = None,
    flags: tuple[str, ...] = DEFAULT_FLAGS,
) -> CompiledFortran:
    source = emit_fortran(nir, backend=backend)
    driver = emit_driver(nir)
    artifact = compile_fortran_source(source, driver=driver, workdir=workdir, flags=flags)
    if artifact.nir_digest != native_ir_digest(nir):
        raise RuntimeError("compiled digest mismatch")
    return artifact


def compile_primitive_driver(
    *,
    workdir: Path | str | None = None,
    runtime: Path | str | None = None,
    flags: tuple[str, ...] = DEFAULT_FLAGS,
) -> Path:
    compiler = find_gfortran()
    runtime_path = Path(runtime) if runtime is not None else RUNTIME_FILE
    if workdir is None:
        workdir = Path(tempfile.mkdtemp(prefix="fepdmc_toolkit_prim_"))
    else:
        workdir = Path(workdir)
        workdir.mkdir(parents=True, exist_ok=True)
    rt_copy = workdir / "native_kernel_runtime_v1.f90"
    drv_copy = workdir / "nk_prim_driver.f90"
    shutil.copyfile(runtime_path, rt_copy)
    shutil.copyfile(PRIM_DRIVER, drv_copy)
    exe = workdir / "nk_prim"
    cmd = [compiler, *flags, "-J", str(workdir), str(rt_copy), str(drv_copy), "-o", str(exe)]
    proc = subprocess.run(cmd, capture_output=True, text=True, check=False)
    if proc.returncode != 0:
        raise RuntimeError(f"primitive driver compile failed: {proc.stderr}")
    return exe


def run_primitive(exe: Path, stdin: str) -> str:
    proc = subprocess.run(
        [str(exe)],
        input=stdin,
        text=True,
        capture_output=True,
        check=False,
    )
    if proc.returncode != 0:
        raise RuntimeError(f"primitive driver failed: {proc.stderr}")
    return proc.stdout
