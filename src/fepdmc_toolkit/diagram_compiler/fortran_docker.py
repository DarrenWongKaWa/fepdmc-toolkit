"""Compile/run generated kernels in the public FEP-DMC gfortran environment."""

from __future__ import annotations

import os
import shutil
import subprocess
from pathlib import Path

from .fortran_compile import DEFAULT_FLAGS, FortranResult

DEFAULT_IMAGE = os.environ.get("FEPDMC_TOOLKIT_DC_DOCKER", "fepdmc-toolkit-u22:local")
PUBLIC_IMAGE = "public.ecr.aws/docker/library/ubuntu:22.04"
FFLAGS = DEFAULT_FLAGS + ("-J", ".")


def docker_available(image: str = DEFAULT_IMAGE) -> bool:
    if shutil.which("docker") is None:
        return False
    probe = subprocess.run(
        ["docker", "image", "inspect", image],
        capture_output=True,
        text=True,
        check=False,
    )
    return probe.returncode == 0


def docker_run(
    image: str, args: list[str], *, workdir: Path, stdin: str | None = None
) -> subprocess.CompletedProcess[str]:
    cmd = ["docker", "run", "--rm"]
    if stdin is not None:
        cmd.append("-i")
    cmd.extend(
        [
            "-v",
            f"{workdir.resolve()}:/work",
            "-w",
            "/work",
            image,
            *args,
        ]
    )
    return subprocess.run(cmd, input=stdin, text=True, capture_output=True, check=False)


def gfortran_version(image: str = DEFAULT_IMAGE) -> str:
    proc = docker_run(image, ["gfortran", "--version"], workdir=Path("/tmp"))
    if proc.returncode != 0:
        raise RuntimeError(proc.stderr)
    return proc.stdout.splitlines()[0]


def compile_in_docker(
    *,
    runtime: Path,
    generated: Path,
    driver: Path,
    workdir: Path,
    image: str = DEFAULT_IMAGE,
    exe_name: str = "dc_validate.x",
) -> Path:
    workdir.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(runtime, workdir / "native_kernel_runtime_v1.f90")
    shutil.copyfile(generated, workdir / "generated_da_kernel.f90")
    shutil.copyfile(driver, workdir / "dc_validate.f90")
    proc = docker_run(
        image,
        [
            "gfortran",
            *FFLAGS,
            "native_kernel_runtime_v1.f90",
            "generated_da_kernel.f90",
            "dc_validate.f90",
            "-o",
            exe_name,
        ],
        workdir=workdir,
    )
    if proc.returncode != 0:
        raise RuntimeError(f"docker gfortran failed: {proc.stderr}\n{proc.stdout}")
    exe = workdir / exe_name
    if not exe.exists():
        raise RuntimeError("docker compile produced no executable")
    return exe


def run_in_docker(exe_name: str, stdin: str, *, workdir: Path, image: str = DEFAULT_IMAGE) -> FortranResult:
    proc = docker_run(image, [f"./{exe_name}"], workdir=workdir, stdin=stdin)
    if proc.returncode != 0:
        raise RuntimeError(f"native kernel exited {proc.returncode}: {proc.stderr}\n{proc.stdout}")
    line = proc.stdout.strip().splitlines()[-1]
    parts = line.split()
    if len(parts) != 12:
        raise RuntimeError(f"unexpected driver output {proc.stdout!r}")
    return FortranResult(
        status=int(parts[0]),
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
