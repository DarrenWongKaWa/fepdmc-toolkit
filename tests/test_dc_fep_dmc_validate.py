"""Stage 6: FEP-DMC-toolchain validation of the stage-5 generated kernel."""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path

import pytest

import fepdmc_toolkit.diagram_compiler.fortran_docker as docker_runtime
from fepdmc_toolkit.diagram_compiler import (
    build_twoband_ir,
    compile_delayed_acceptance,
    emit_fortran,
    lower_native_kernel,
)
from fepdmc_toolkit.diagram_compiler.fortran_compile import ATOL, RTOL, RUNTIME_FILE
from fepdmc_toolkit.diagram_compiler.fortran_docker import (
    DEFAULT_IMAGE,
    compile_in_docker,
    docker_available,
    run_in_docker,
)
from fepdmc_toolkit.diagram_compiler.fortran_fixtures import twoband_oracle_cases

ROOT = Path(__file__).resolve().parents[1]
GOLDEN = json.loads((ROOT / "validation/diagram_compiler/golden.json").read_text())
PIN = json.loads((ROOT / "provenance/UPSTREAM_FEP_DMC.json").read_text())
FEP_DEFAULT = Path(os.environ["FEPDMC_TOOLKIT_FEP_DMC"]) if os.environ.get("FEPDMC_TOOLKIT_FEP_DMC") else None


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def test_docker_availability_requires_the_configured_local_image(monkeypatch):
    """Docker installed without the pinned image must not run native tests."""

    class Result:
        def __init__(self, returncode):
            self.returncode = returncode

    seen = []

    def inspect(cmd, **kwargs):
        seen.append(cmd)
        return Result(1)

    monkeypatch.setattr(docker_runtime.shutil, "which", lambda name: "/usr/bin/docker")
    monkeypatch.setattr(docker_runtime.subprocess, "run", inspect)
    assert docker_runtime.docker_available("fepdmc-toolkit-u22:missing") is False
    assert seen == [["docker", "image", "inspect", "fepdmc-toolkit-u22:missing"]]


def test_generated_source_and_runtime_match_the_validated_artifacts():
    kernel = compile_delayed_acceptance(build_twoband_ir(1))
    nir = lower_native_kernel(kernel)
    src = emit_fortran(nir)
    assert hashlib.sha256(src.encode("utf-8")).hexdigest() == GOLDEN["fortran_source_sha256"]
    assert _sha(RUNTIME_FILE) == GOLDEN["runtime_sha256"]
    artifact = ROOT / "validation/diagram_compiler/generated_twoband_n1.f90"
    assert _sha(artifact) == GOLDEN["fortran_source_sha256"]


def test_backend_still_does_not_import_diagramir():
    import ast

    tree = ast.parse((ROOT / "src/fepdmc_toolkit/diagram_compiler/fortran_codegen.py").read_text())
    imported = []
    for node in tree.body:
        if isinstance(node, ast.ImportFrom):
            imported.append(node.module)
    assert "ir" not in imported
    assert "da_kernel" not in imported
    assert "cheap_policy" not in imported


def test_fep_dmc_pin_pristine():
    import sys

    sys.path.insert(0, str(ROOT / "integration" / "diagram_compiler"))
    import dc_lib as dc

    if FEP_DEFAULT is None or not FEP_DEFAULT.is_dir():
        pytest.skip("pinned FEP-DMC clone not present")
    pin = dc.load_json(dc.DEFAULT_PIN)
    state, detail = dc.classify(FEP_DEFAULT, pin)
    assert state == dc.STATE_PUBLIC_PRISTINE, detail
    assert dc.compiler_kernel_needles_in_tree(FEP_DEFAULT) == []
    assert dc.git_head(FEP_DEFAULT) == PIN["reference_commit"]


@pytest.mark.skipif(not docker_available(), reason="configured Docker image not available")
def test_ubuntu_gfortran_kernel_matches_oracle(tmp_path):
    from fepdmc_toolkit.diagram_compiler import emit_driver

    kernel = compile_delayed_acceptance(build_twoband_ir(1))
    nir = lower_native_kernel(kernel)
    generated = tmp_path / "generated_da_kernel.f90"
    driver = tmp_path / "dc_validate.f90"
    generated.write_text(emit_fortran(nir), newline="\n")
    driver.write_text(emit_driver(nir), newline="\n")
    workdir = tmp_path / "build"
    compile_in_docker(runtime=RUNTIME_FILE, generated=generated, driver=driver, workdir=workdir, image=DEFAULT_IMAGE)
    cases = twoband_oracle_cases(nir, kernel)
    for case in cases:
        ft = run_in_docker("dc_validate.x", case.text, workdir=workdir, image=DEFAULT_IMAGE)
        exp = case.expected
        assert ft.status == exp["status"]
        if exp["status"] != 0:
            assert ft.accepted is False
            continue
        assert ft.accepted == exp["accepted"]
        assert ft.reject_stage == exp["reject_stage"]
        assert ft.exact_y_evaluated == exp["exact_y_evaluated"]
        if exp["ell_hat"] is not None:
            assert ft.ell_hat == pytest.approx(exp["ell_hat"], abs=ATOL, rel=RTOL)
        if exp["ell_R"] is None:
            assert ft.ell_R_valid is False
        else:
            assert ft.ell_R_valid is True
            assert ft.ell_R == pytest.approx(exp["ell_R"], abs=ATOL, rel=RTOL)
        if case.name == "stage1_reject":
            assert ft.n_exact_graph == 0
            assert ft.n_electron_twoband == 0
            assert ft.n_trace == 0
        if case.name in ("accept", "stage2_reject", "asymmetric_q", "uncached_exact_x"):
            assert ft.n_exact_graph > 0
        if case.name == "cached_exact_x":
            assert ft.n_exact_graph == 1
