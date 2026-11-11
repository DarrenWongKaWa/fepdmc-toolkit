"""Stage 4: NativeKernelIR, the typed control-flow form of the delayed-acceptance kernel.

The two-band n = 1 kernel is lowered to NativeKernelIR. The IR must validate, must
equal the recorded IR (validation/diagram_compiler/NATIVE_KERNEL_IR.json), must not
reach an exact operation on the stage-1 reject path, and its interpreter must agree
with the Python kernel on every path and keep detailed balance on the three-state
chain.

Records: NATIVE_KERNEL_IR.json, binding_layout.json, primitive_inventory.json,
cfg_paths.json, differential_results.csv, finite_state_native.csv and summary.json.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np
from _common import (
    BALANCE_TOL,
    CHEAP,
    PI,
    Q_ASYMMETRIC,
    RECORDS,
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
    interpret_native_kernel,
    lower_native_kernel,
    native_acceptance_probability,
)
from fepdmc_toolkit.diagram_compiler.fortran_codegen import native_ir_digest
from fepdmc_toolkit.diagram_compiler.native_validate import (
    exact_candidate_ops_on_reject_path,
    reachable,
    stage1_reject_blocks,
    validate_native_ir,
)


def main(output: Path) -> int:
    output.mkdir(parents=True, exist_ok=True)
    kernel = compile_delayed_acceptance(build_twoband_ir(1))
    nir = lower_native_kernel(kernel)
    validate_native_ir(nir)
    ir_json = json.loads(nir.to_json())
    write_json(output / "NATIVE_KERNEL_IR.json", ir_json)
    write_json(output / "binding_layout.json", nir.layout.to_dict())
    write_json(
        output / "primitive_inventory.json",
        {
            "cheap_ops": sorted({op.op for op in nir.graph("cheap").ops}),
            "control_ops": sorted({op.op for block in nir.blocks for op in block.ops}),
            "exact_ops": sorted({op.op for op in nir.graph("exact").ops}),
            "primitive_set": "native_kernel_v1",
            "schema_version": 1,
            "signatures": {
                "prim.electron_scalar": "(k_eff:f64, t:f64, dtau:f64) -> c128; Gel=exp(-2t(1-cos k_eff) dtau)",
                "prim.electron_twoband": "(k_eff,t,dtau,delta,gap) -> c128_m2; U=exp(-H dtau) via eigh",
                "prim.phonon": "(omega:f64, dtau:f64) -> c128; D=exp(-omega dtau)",
                "prim.prefactor": "(g:f64, L:f64; n:i64 attr) -> c128; (g/sqrt(L))**(2n)",
                "prim.vertex_sigmaz": "(g:f64, L:f64) -> c128_m2; (g/sqrt(L))*diag(1,-1)",
            },
        },
    )
    exact_on_reject = exact_candidate_ops_on_reject_path(nir)
    write_json(
        output / "cfg_paths.json",
        {
            "blocks": [block.id for block in nir.blocks],
            "entry": nir.entry,
            "exact_ops_on_stage1_reject": exact_on_reject,
            "reachable_from_entry": sorted(reachable(nir, nir.entry)),
            "stage1_reject_blocks": sorted(stage1_reject_blocks(nir)),
        },
    )
    x, y = twoband_binding("x"), twoband_binding("y")
    py_probe = kernel.evaluate_transition(x, y, 0.0, 1e-16, 1e-16)
    alpha2 = math.exp(min(0.0, py_probe.ell_R - py_probe.ell_hat))
    cases = [
        ("stage1_reject", 0.0, 0.9, 0.1),
        ("stage2_reject", 0.0, 1e-16, (alpha2 + 1.0) / 2.0),
        ("accept", 0.0, 1e-16, 1e-16),
        ("asymmetric_q", 0.4, 1e-16, 1e-16),
    ]
    rows = []
    for name, log_q, u1, u2 in cases:
        py = kernel.evaluate_transition(x, y, log_q, u1, u2)
        native = interpret_native_kernel(nir, x=x, y=y, log_q=log_q, u1=u1, u2=u2)
        rows.append(
            {
                "accepted_match": py.accepted == native.accepted,
                "ell_R_native": native.ell_R,
                "ell_hat_native": native.ell_hat,
                "exact_y_match": py.exact_y_evaluated == native.exact_y_evaluated,
                "name": name,
                "native_accepted": native.accepted,
                "py_accepted": py.accepted,
            }
        )
    write_csv(output / "differential_results.csv", rows)
    n = len(PI)
    p = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            log_q = math.log(Q_ASYMMETRIC[j, i]) - math.log(Q_ASYMMETRIC[i, j])
            a = native_acceptance_probability(
                math.log(CHEAP[i]),
                math.log(CHEAP[j]),
                math.log(PI[i]),
                math.log(PI[j]),
                log_q,
            )
            p[i, j] = Q_ASYMMETRIC[i, j] * a
        p[i, i] = 1.0 - p[i].sum()
    residual, stat = balance_residuals(p)
    write_csv(
        output / "finite_state_native.csv",
        [
            {
                "max_balance_residual": residual,
                "max_stationarity_residual": stat,
                "status": "PASS" if residual <= BALANCE_TOL else "FAIL",
            }
        ],
    )
    return finish(
        "stage 4 (NativeKernelIR)",
        output,
        {
            "ir_equals_the_recorded_ir": ir_json == json.loads((RECORDS / "NATIVE_KERNEL_IR.json").read_text())
            and native_ir_digest(nir) == golden()["native_ir_digest"],
            "no_exact_operation_on_the_stage1_reject_path": exact_on_reject == [],
            "interpreter_agrees_with_the_python_kernel": all(r["accepted_match"] and r["exact_y_match"] for r in rows),
            "detailed_balance_on_the_three_state_chain": residual <= BALANCE_TOL,
        },
        max_balance_residual=residual,
        native_ir_digest=native_ir_digest(nir),
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--output", type=Path, default=ROOT / "build/validation/stage4_native_ir")
    sys.exit(main(parser.parse_args().output))
