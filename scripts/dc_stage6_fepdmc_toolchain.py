"""Stage 6: the generated kernel in the Ubuntu 22.04 gfortran image, next to the pinned FEP-DMC tree.

The emitted two-band n = 1 kernel, the primitive runtime and an emitted driver are
compiled with gfortran inside the Docker image built from
integration/diagram_compiler/Dockerfile.u22 and run on fixtures made by the
NativeKernelIR interpreter: stage-1 reject, stage-2 reject, accept, asymmetric
proposal, cached and uncached exact(x), and the error statuses. Discrete fields
must match exactly, the scores within the tolerance. The compiled score kernel must
keep detailed balance on the three-state chain.

The pinned FEP-DMC clone is a witness only: it must be pristine and must not
contain the compiler kernel; it is not modified, and production FEP-DMC does not
call the generated kernel.

Records: generated_twoband_n1.f90, dc_validate.f90, Makefile, fep_dmc_pin.json,
fixtures/, fixture_manifest.json, native_differential.csv, native_call_counts.csv,
finite_state_native.csv, build_environment.json, integration_provenance.json and
summary.json.
"""

from __future__ import annotations

import argparse
import hashlib
import math
import os
import sys
import tempfile
from pathlib import Path

import numpy as np
from _common import BALANCE_TOL, CHEAP, PI, Q_ASYMMETRIC, ROOT, balance_residuals, finish, golden, write_csv, write_json

from fepdmc_toolkit.diagram_compiler import (
    build_twoband_ir,
    compile_delayed_acceptance,
    emit_driver,
    emit_fortran,
    interpret_native_kernel,
    lower_native_kernel,
    lower_native_score_kernel,
)
from fepdmc_toolkit.diagram_compiler.fortran_codegen import native_ir_digest
from fepdmc_toolkit.diagram_compiler.fortran_compile import ATOL, RTOL, RUNTIME_FILE
from fepdmc_toolkit.diagram_compiler.fortran_docker import (
    DEFAULT_IMAGE,
    FFLAGS,
    PUBLIC_IMAGE,
    compile_in_docker,
    docker_available,
    docker_run,
    gfortran_version,
    run_in_docker,
)
from fepdmc_toolkit.diagram_compiler.fortran_fixtures import (
    binding_to_dict,
    twoband_oracle_cases,
)

sys.path.insert(0, str(ROOT / "integration" / "diagram_compiler"))
import dc_lib as dc  # noqa: E402


def sha256_text(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def sha256_path(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def image_platform(image: str) -> str:
    """Operating system and machine of the image, for example 'Ubuntu 22.04.5 LTS x86_64'."""
    proc = docker_run(image, ["sh", "-c", '. /etc/os-release && echo "$PRETTY_NAME $(uname -m)"'], workdir=Path("/tmp"))
    return proc.stdout.strip() if proc.returncode == 0 and proc.stdout.strip() else "unknown"


def match_row(name: str, expected: dict, ft) -> dict:
    ell_hat_err = None
    ell_r_err = None
    if expected["status"] == 0 and expected["ell_hat"] is not None:
        ell_hat_err = abs(ft.ell_hat - expected["ell_hat"])
    if expected["ell_R"] is not None:
        ell_r_err = abs(ft.ell_R - expected["ell_R"])
    return {
        "name": name,
        "status_match": expected["status"] == ft.status,
        "accepted_match": expected["accepted"] == ft.accepted,
        "reject_stage_match": expected["reject_stage"] == ft.reject_stage,
        "exact_y_match": expected["exact_y_evaluated"] == ft.exact_y_evaluated,
        "ell_hat_abs_err": ell_hat_err,
        "ell_R_abs_err": ell_r_err,
        "py_status": expected["status"],
        "ft_status": ft.status,
        "py_accepted": expected["accepted"],
        "ft_accepted": ft.accepted,
        "py_reject_stage": expected["reject_stage"],
        "ft_reject_stage": ft.reject_stage,
        "py_ell_hat": expected["ell_hat"],
        "ft_ell_hat": ft.ell_hat if expected["status"] == 0 else None,
        "py_ell_R": expected["ell_R"],
        "ft_ell_R": ft.ell_R if ft.ell_R_valid else None,
        "n_exact_graph": ft.n_exact_graph,
        "n_electron_twoband": ft.n_electron_twoband,
        "n_trace": ft.n_trace,
    }


def close(value: float, reference: float) -> bool:
    return abs(value - reference) <= ATOL + RTOL * abs(reference)


def agrees(expected: dict, ft) -> bool:
    """The native result equals the interpreter's: discrete fields exactly, scores within the tolerance."""
    if ft.status != expected["status"]:
        return False
    if expected["status"] != 0:
        return ft.accepted is False
    if (
        ft.accepted != expected["accepted"]
        or ft.reject_stage != expected["reject_stage"]
        or ft.exact_y_evaluated != expected["exact_y_evaluated"]
    ):
        return False
    if expected["ell_hat"] is not None and not close(ft.ell_hat, expected["ell_hat"]):
        return False
    if expected["ell_R"] is None:
        return ft.ell_R_valid is False
    return ft.ell_R_valid is True and close(ft.ell_R, expected["ell_R"])


def main(output: Path, fep_dmc: Path, image: str) -> int:
    output.mkdir(parents=True, exist_ok=True)
    reference = golden()
    kernel = compile_delayed_acceptance(build_twoband_ir(1))
    nir = lower_native_kernel(kernel)
    generated = emit_fortran(nir)
    driver = emit_driver(nir)
    gen_sha = sha256_text(generated)
    artifacts_ok = (
        gen_sha == reference["fortran_source_sha256"]
        and native_ir_digest(nir) == reference["native_ir_digest"]
        and sha256_path(RUNTIME_FILE) == reference["runtime_sha256"]
    )
    (output / "generated_twoband_n1.f90").write_text(generated, encoding="utf-8", newline="\n")
    (output / "dc_validate.f90").write_text(driver, encoding="utf-8", newline="\n")
    makefile = (ROOT / "integration/diagram_compiler/Makefile").read_text(encoding="utf-8")
    (output / "Makefile").write_text(makefile, encoding="utf-8", newline="\n")

    pin = dc.load_json(dc.DEFAULT_PIN)
    state, detail = dc.classify(fep_dmc, pin)
    needles = dc.compiler_kernel_needles_in_tree(fep_dmc)
    write_json(
        output / "fep_dmc_pin.json",
        {
            "detail": detail,
            "head": dc.git_head(fep_dmc),
            "needles": needles,
            "production_hashes": dc.production_hashes(fep_dmc, pin),
            "reference_commit": pin["reference_commit"],
            "state": state,
        },
    )
    pin_ok = state == dc.STATE_PUBLIC_PRISTINE and not needles
    if not artifacts_ok or not pin_ok:
        # Nothing is compiled against artifacts or a tree that are not the validated ones.
        return finish(
            "stage 6 (FEP-DMC toolchain)",
            output,
            {
                "emitted_artifacts_equal_the_validated_ones": artifacts_ok,
                "fep_dmc_tree_is_pristine_and_free_of_the_kernel": pin_ok,
            },
            fep_dmc_state=f"{state}: {detail}",
        )

    cases = twoband_oracle_cases(nir, kernel)
    manifest = []
    text_dir = output / "fixtures"
    text_dir.mkdir(exist_ok=True)
    for case in cases:
        path = text_dir / f"{case.name}.txt"
        path.write_text(case.text, encoding="utf-8", newline="\n")
        expected = dict(case.expected)
        for key in ("ell_hat", "ell_R"):
            value = expected[key]
            if isinstance(value, float) and not math.isfinite(value):
                expected[key] = str(value)
        log_q = case.log_q if math.isfinite(case.log_q) else str(case.log_q)
        manifest.append(
            {
                "expected": expected,
                "exact_x_valid": case.exact_x_valid,
                "exact_log_weight_x": case.exact_log_weight_x,
                "log_q": log_q,
                "name": case.name,
                "text_sha256": sha256_path(path),
                "u1": case.u1,
                "u2": case.u2,
                "x": binding_to_dict(case.x, nir.layout),
                "y": binding_to_dict(case.y, nir.layout),
            }
        )
    write_json(output / "fixture_manifest.json", {"cases": manifest, "schema": "diagram_compiler_stage6_fixtures_v1"})

    build = Path(tempfile.mkdtemp(prefix="stage6_u22_"))
    compile_in_docker(
        runtime=RUNTIME_FILE,
        generated=output / "generated_twoband_n1.f90",
        driver=output / "dc_validate.f90",
        workdir=build,
        image=image,
    )
    version = gfortran_version(image)
    rows = []
    call_rows = []
    fixtures_ok = True
    for case in cases:
        ft = run_in_docker("dc_validate.x", case.text, workdir=build, image=image)
        fixtures_ok = fixtures_ok and agrees(case.expected, ft)
        rows.append(match_row(case.name, case.expected, ft))
        call_rows.append(
            {
                "name": case.name,
                "n_electron_twoband": ft.n_electron_twoband,
                "n_exact_graph": ft.n_exact_graph,
                "n_matmul": ft.n_matmul,
                "n_trace": ft.n_trace,
                "n_vertex_sigmaz": ft.n_vertex_sigmaz,
                "stage1_reject": case.name == "stage1_reject",
            }
        )
    write_csv(output / "native_differential.csv", rows)
    write_csv(output / "native_call_counts.csv", call_rows)

    score = lower_native_score_kernel()
    score_src = emit_fortran(score)
    score_drv = emit_driver(score)
    score_dir = Path(tempfile.mkdtemp(prefix="stage6_score_"))
    (score_dir / "generated_score.f90").write_text(score_src, newline="\n")
    (score_dir / "score_driver.f90").write_text(score_drv, newline="\n")
    compile_in_docker(
        runtime=RUNTIME_FILE,
        generated=score_dir / "generated_score.f90",
        driver=score_dir / "score_driver.f90",
        workdir=score_dir,
        image=image,
        exe_name="dc_score.x",
    )
    n = len(PI)
    p = np.zeros((n, n))
    score_status_ok = True
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            log_q = math.log(Q_ASYMMETRIC[j, i]) - math.log(Q_ASYMMETRIC[i, j])
            stdin = (
                f"{math.log(CHEAP[i])} {math.log(CHEAP[j])} {math.log(PI[i])} {math.log(PI[j])}\n"
                f"{log_q} {1e-16} {1e-16}\n"
            )
            ft = run_in_docker("dc_score.x", stdin, workdir=score_dir, image=image)
            py = interpret_native_kernel(
                score,
                log_wx=math.log(CHEAP[i]),
                log_wy=math.log(CHEAP[j]),
                log_px=math.log(PI[i]),
                log_py=math.log(PI[j]),
                log_q=log_q,
                u1=1e-16,
                u2=1e-16,
            )
            score_status_ok = score_status_ok and ft.status == 0 and py.status == 0
            a1 = 1.0 if ft.ell_hat >= 0.0 else math.exp(ft.ell_hat)
            delta = ft.ell_R - ft.ell_hat
            a2 = 1.0 if delta >= 0.0 else math.exp(delta)
            p[i, j] = Q_ASYMMETRIC[i, j] * a1 * a2
        p[i, i] = 1.0 - p[i].sum()
    residual, stat = balance_residuals(p)
    write_csv(
        output / "finite_state_native.csv",
        [{"balance_residual": residual, "stationarity_residual": stat, "atol": ATOL, "rtol": RTOL}],
    )

    exe = build / "dc_validate.x"
    write_json(
        output / "build_environment.json",
        {
            "docker_image": image,
            "fflags": list(FFLAGS),
            "gfortran_version": version,
            "os": image_platform(image),
            "base_image": PUBLIC_IMAGE,
            "qe_hdf5_mpi_linked": False,
        },
    )
    write_json(
        output / "integration_provenance.json",
        {
            "validation_executable": "separate from the FEP-DMC tree",
            "binary_sha256": sha256_path(exe),
            "docker_image": image,
            "driver_sha256": sha256_text(driver),
            "fep_dmc_commit": pin["reference_commit"],
            "fortran_source_sha256": gen_sha,
            "gfortran_version": version,
            "native_ir_digest": native_ir_digest(nir),
            "production_fep_dmc_modified": False,
            "runtime_sha256": sha256_path(RUNTIME_FILE),
            "cheap_policy": "propagator_only_v1",
            "primitive_set": "native_kernel_v1",
            "target": "positive_real_F_v1",
            "tolerances": {"atol": ATOL, "discrete": "exact", "rtol": RTOL},
        },
    )
    stage1 = next(r for r in call_rows if r["stage1_reject"])
    errors = [
        e
        for r in rows
        if r["py_status"] == 0
        for e in (r["ell_hat_abs_err"], r["ell_R_abs_err"])
        if e is not None and math.isfinite(e)
    ]
    return finish(
        "stage 6 (FEP-DMC toolchain)",
        output,
        {
            "emitted_artifacts_equal_the_validated_ones": artifacts_ok,
            "fep_dmc_tree_is_pristine_and_free_of_the_kernel": pin_ok,
            "native_kernel_matches_the_interpreter_on_every_fixture": fixtures_ok,
            "no_exact_primitive_on_the_stage1_reject_path": stage1["n_exact_graph"] == 0
            and stage1["n_electron_twoband"] == 0
            and stage1["n_trace"] == 0,
            "score_kernel_returns_ok": score_status_ok,
            "detailed_balance_on_the_three_state_chain": residual <= BALANCE_TOL and stat <= BALANCE_TOL,
        },
        fixtures=len(rows),
        gfortran_version=version,
        max_balance_residual=residual,
        max_score_abs_err=max(errors) if errors else 0.0,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--output", type=Path, default=ROOT / "build/validation/stage6_fepdmc_toolchain")
    parser.add_argument(
        "--fep-dmc",
        type=Path,
        default=(Path(os.environ["FEPDMC_TOOLKIT_FEP_DMC"]) if os.environ.get("FEPDMC_TOOLKIT_FEP_DMC") else None),
        help="clone of the pinned FEP-DMC commit (default: FEPDMC_TOOLKIT_FEP_DMC)",
    )
    parser.add_argument("--image", default=DEFAULT_IMAGE, help=f"Docker image (default: {DEFAULT_IMAGE})")
    args = parser.parse_args()
    if args.fep_dmc is None:
        parser.error("--fep-dmc or FEPDMC_TOOLKIT_FEP_DMC is required; there is no default path")
    if not docker_available(args.image):
        parser.error(
            f"Docker image {args.image!r} not found; build it with "
            "docker build -f integration/diagram_compiler/Dockerfile.u22 -t fepdmc-toolkit-u22:local ."
        )
    sys.exit(main(args.output, args.fep_dmc, args.image))
