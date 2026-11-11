"""Stage 3: the two-stage delayed-acceptance kernel (Python).

The two-band n = 1 kernel is driven through its three paths (stage-1 reject,
stage-2 reject, accept), and the kernel algebra is checked on a three-state chain
with a symmetric and an asymmetric proposal: detailed balance and stationarity for
the exact target must hold to 1e-12 although the cheap weight differs from it.

Records: DA_KERNEL_SPEC.json, TARGET_POLICY.json, PROPOSAL_SPEC.json,
demo_transitions.json, call_counts.csv, finite_state_balance.csv and summary.json.
"""

from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

from _common import (
    BALANCE_TOL,
    CHEAP,
    PI,
    Q_ASYMMETRIC,
    Q_SYMMETRIC,
    ROOT,
    balance_residuals,
    finish,
    twoband_binding,
    write_csv,
    write_json,
)

from fepdmc_toolkit.diagram_compiler import (
    Binding,
    build_twoband_ir,
    compile_delayed_acceptance,
    positive_real_F_v1,
    transition_matrix,
)
from fepdmc_toolkit.diagram_compiler.proposal import ProposalSpec


def result_row(name: str, x: Binding, y: Binding, log_q: float, u1: float, u2: float, result) -> dict:
    return {
        "accepted": result.accepted,
        "cheap_log_weight_x": result.cheap_log_weight_x,
        "cheap_log_weight_y": result.cheap_log_weight_y,
        "ell_R": result.ell_R,
        "ell_hat": result.ell_hat,
        "exact_log_weight_x": result.exact_log_weight_x,
        "exact_log_weight_y": result.exact_log_weight_y,
        "exact_y_evaluated": result.exact_y_evaluated,
        "g_x": x.g,
        "g_y": y.g,
        "k_x": x.k,
        "k_y": y.k,
        "log_q_reverse_minus_forward": log_q,
        "name": name,
        "stage1_pass": result.stage1_pass,
        "stage2_pass": result.stage2_pass,
        "u1": u1,
        "u2": u2,
    }


def main(output: Path) -> int:
    output.mkdir(parents=True, exist_ok=True)
    ir = build_twoband_ir(1)
    kernel = compile_delayed_acceptance(ir, proposal_spec=ProposalSpec.provided_log_ratio())
    write_json(output / "DA_KERNEL_SPEC.json", kernel.spec.to_dict())
    write_json(output / "TARGET_POLICY.json", positive_real_F_v1().to_dict())
    write_json(output / "PROPOSAL_SPEC.json", ProposalSpec.provided_log_ratio().to_dict())
    x, y = twoband_binding("x"), twoband_binding("y")
    probe = kernel.evaluate_transition(x, y, 0.0, 1e-16, 1e-16)
    alpha2 = 1.0 if probe.ell_R - probe.ell_hat >= 0 else math.exp(probe.ell_R - probe.ell_hat)
    u2_reject = (alpha2 + 1.0) / 2.0
    cases = [
        ("stage1_reject", x, y, 0.0, 0.9, 0.1),
        ("stage2_reject", x, y, 0.0, 1e-16, u2_reject),
        ("accept", x, y, 0.0, 1e-16, 1e-16),
    ]
    rows = []
    n_s1_reject = n_s1_pass = n_exact_y = n_s2_reject = n_accept = 0
    for name, cx, cy, log_q, u1, u2 in cases:
        result = kernel.evaluate_transition(cx, cy, log_q, u1, u2)
        rows.append(result_row(name, cx, cy, log_q, u1, u2, result))
        if not result.stage1_pass:
            n_s1_reject += 1
        else:
            n_s1_pass += 1
            n_exact_y += int(result.exact_y_evaluated)
            if result.stage2_pass:
                n_accept += 1
            else:
                n_s2_reject += 1
    write_json(output / "demo_transitions.json", {"family": ir.family.family_id, "transitions": rows})
    write_csv(
        output / "call_counts.csv",
        [
            {
                "accepts": n_accept,
                "exact_candidate_evaluations": n_exact_y,
                "stage1_passes": n_s1_pass,
                "stage1_rejects": n_s1_reject,
                "stage2_rejects": n_s2_reject,
                "transitions_attempted": len(cases),
            }
        ],
    )
    balance_rows = []
    for name, q in (("symmetric", Q_SYMMETRIC), ("asymmetric", Q_ASYMMETRIC)):
        residual, stat = balance_residuals(transition_matrix(PI, CHEAP, q))
        balance_rows.append(
            {
                "max_balance_residual": residual,
                "max_stationarity_residual": stat,
                "proposal": name,
                "status": "PASS" if residual <= BALANCE_TOL and stat <= BALANCE_TOL else "FAIL",
            }
        )
    write_csv(output / "finite_state_balance.csv", balance_rows)
    by_name = {row["name"]: row for row in rows}
    return finish(
        "stage 3 (delayed-acceptance kernel)",
        output,
        {
            "stage1_reject_skips_the_exact_candidate": not by_name["stage1_reject"]["stage1_pass"]
            and not by_name["stage1_reject"]["exact_y_evaluated"],
            "stage2_reject_path": by_name["stage2_reject"]["stage1_pass"]
            and not by_name["stage2_reject"]["stage2_pass"]
            and not by_name["stage2_reject"]["accepted"],
            "accept_path": bool(by_name["accept"]["accepted"]),
            "exact_candidate_evaluated_once_per_stage1_pass": n_exact_y == n_s1_pass,
            "detailed_balance_on_the_three_state_chain": all(r["status"] == "PASS" for r in balance_rows),
        },
        demo_paths=[row["name"] for row in rows],
        max_balance_residual=max(r["max_balance_residual"] for r in balance_rows),
        target_policy="positive_real_F_v1",
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--output", type=Path, default=ROOT / "build/validation/stage3_da_kernel")
    sys.exit(main(parser.parse_args().output))
