"""Stage 1: the exact compiler against the R1 reference and an independent index oracle.

For scalar Holstein (n = 1-4) and Hermitian two-band (n = 1-3), sixteen bindings on
each of two base points are evaluated diagram by diagram, through the plain and the
shared DAG, through the grouped lowering and through the compiled evaluator. Every
result must agree with the R1 recursive builder and with the explicit enumeration
within ``numeric_close`` (atol 1e-12, rtol 1e-10).

Records: differential.csv, graph_stats.csv, R1_ORACLE.json, a worked example
(r1_scalar_n2.json, binding.json, demo_result.json) and summary.json.
"""

from __future__ import annotations

import argparse
import hashlib
import sys
from collections import Counter
from dataclasses import asdict, replace
from pathlib import Path

import numpy as np
from _common import ROOT, finish, golden, write_csv, write_json

from fepdmc_toolkit.diagram_compiler import (
    DiagramIR,
    compile_evaluator,
    evaluate_dag,
    evaluate_diagramwise,
    import_r1_binding,
    import_r1_spec,
    lower_diagramwise,
    lower_grouped,
    numeric_close,
)
from fepdmc_toolkit.scheme_d.exec import interpret
from fepdmc_toolkit.scheme_d.physics_first.contract import HolsteinRing, TwoBandHolstein
from fepdmc_toolkit.scheme_d.physics_first.index_oracle import group_sum_matrix_twoband, group_sum_scalar
from fepdmc_toolkit.scheme_d.physics_first.symbolic.bindings import apply_scenario, base_binding
from fepdmc_toolkit.scheme_d.physics_first.symbolic.plan_spec import scalar_spec, twoband_spec
from fepdmc_toolkit.scheme_d.physics_first.symbolic.r1 import r1_symbolic_build
from fepdmc_toolkit.scheme_d.physics_first.symbolic.r1_twoband import r1_twoband_build, structural_eye
from fepdmc_toolkit.scheme_d.physics_first.symbolic.reference_binder import ReferenceBinder

# Sources of the two reference implementations; their hashes go into R1_ORACLE.json.
ORACLE_SOURCES = (
    "src/fepdmc_toolkit/scheme_d/physics_first/contract.py",
    "src/fepdmc_toolkit/scheme_d/physics_first/index_oracle.py",
    "src/fepdmc_toolkit/scheme_d/physics_first/fock_witness.py",
    "src/fepdmc_toolkit/scheme_d/physics_first/interface.py",
    "src/fepdmc_toolkit/scheme_d/physics_first/r1_open_dp.py",
    "src/fepdmc_toolkit/scheme_d/physics_first/symbolic/r1.py",
    "src/fepdmc_toolkit/scheme_d/physics_first/symbolic/r1_twoband.py",
    "src/fepdmc_toolkit/scheme_d/physics_first/symbolic/leaves.py",
    "src/fepdmc_toolkit/scheme_d/physics_first/symbolic/reference_binder.py",
    "src/fepdmc_toolkit/scheme_d/physics_first/symbolic/plan_spec.py",
    "src/fepdmc_toolkit/scheme_d/physics_first/symbolic/guards.py",
    "src/fepdmc_toolkit/scheme_d/physics_first/symbolic/bindings.py",
    "src/fepdmc_toolkit/scheme_d/ir.py",
    "src/fepdmc_toolkit/scheme_d/exec.py",
    "src/fepdmc_toolkit/scheme_d/numeric.py",
    "src/fepdmc_toolkit/scheme_d/oracle.py",
)

SCENARIOS = (
    "base",
    "change_k",
    "change_q",
    "change_tau",
    "change_t",
    "change_omega",
    "change_g",
    "g_zero",
    "joint_parameters",
    "joint_sample",
    "repeat_q_then_separate",
    "repeat_dtau_then_separate",
    "change_L_analytic",
    "small_nonzero_g",
    "longer_legal_intervals",
    "same_input_repeat",
)


def main(output: Path) -> int:
    output.mkdir(parents=True, exist_ok=True)
    rows, stats = [], []
    round_trip_ok = True
    for model, orders in [("scalar", range(1, 5)), ("twoband", range(1, 4))]:
        for n in orders:
            spec = scalar_spec(n) if model == "scalar" else twoband_spec(n)
            ir = import_r1_spec(spec)
            revived = DiagramIR.from_json(ir.to_json())
            round_trip_ok = round_trip_ok and revived == ir
            evaluator = compile_evaluator(revived)
            dags = [lower_diagramwise(revived), lower_diagramwise(revived, share=True), lower_grouped(revived)]
            plan = r1_symbolic_build(spec) if model == "scalar" else r1_twoband_build(spec)
            binder = ReferenceBinder(spec, plan.leaves)
            fanout = Counter(i for node in dags[1].nodes for i in node.inputs)
            stats.append(
                {
                    "model": model,
                    "n": n,
                    "diagrams": len(ir.diagrams),
                    "naive_nodes": len(dags[0].nodes),
                    "unique_dag_nodes": len(dags[1].nodes),
                    "nodes_eliminated": len(dags[0].nodes) - len(dags[1].nodes),
                    "shared_nodes": sum(v > 1 for v in fanout.values()),
                    "shared_arithmetic_nodes": sum(
                        v > 1 and dags[1].nodes[i].op not in ("input", "const") for i, v in fanout.items()
                    ),
                    "estimated_operations": sum(node.op not in ("input", "const") for node in dags[1].nodes),
                    "recovered_dp_nodes": len(dags[2].nodes),
                    "r1_nodes": plan.report.dag_nodes,
                    "comparison": "COUNTS_ONLY_NOT_ISOMORPHISM_OR_SPEEDUP",
                }
            )
            for split in ("DEV", "RESERVED_CONFIRMATION"):
                base = base_binding(spec, split)
                for scenario in SCENARIOS:
                    rt = base if scenario == "base" else apply_scenario(base, spec, scenario, 20260918 + n)
                    b = import_r1_binding(rt)
                    env = binder.bind(rt)
                    if model == "twoband":
                        env.update(structural_eye())
                    reference = interpret(plan.graph, env)
                    naive = evaluate_diagramwise(revived, b)
                    outputs = [evaluate_dag(dag, b, revived) for dag in dags] + [evaluator(b)]
                    if model == "scalar":
                        physical = HolsteinRing(L=rt.L, t=rt.t, omega=rt.omega, g=rt.g)
                        independent = {"F": group_sum_scalar(n, rt.tau_array(), rt.q_array(), rt.k, physical)}
                    else:
                        physical = TwoBandHolstein(L=rt.L, t=rt.t, omega=rt.omega, g=rt.g, delta=rt.delta, gap=rt.gap)
                        op = group_sum_matrix_twoband(n, rt.tau_array(), rt.q_array(), rt.k, physical)
                        independent = {"F": np.trace(op), "Op": op}
                    errors = []
                    agrees = True
                    for out in [naive, *outputs]:
                        for key in independent:
                            agrees = agrees and bool(numeric_close(out[key], reference[key]))
                            agrees = agrees and bool(numeric_close(out[key], independent[key]))
                            errors.append(float(np.max(np.abs(out[key] - reference[key]))))
                    rows.append(
                        {
                            "case": f"{spec.spec_id}/{split}/{scenario}",
                            "r1_result": repr(complex(np.asarray(reference["F"]).item())),
                            "naive_generated_result": repr(naive["F"]),
                            "dag_result": repr(outputs[1]["F"]),
                            "max_error_including_Op": max(errors),
                            "status": "PASS" if agrees else "FAIL",
                        }
                    )
            if model == "scalar" and n == 2:
                (output / "r1_scalar_n2.json").write_text(ir.to_json() + "\n")
                rt = replace(
                    base_binding(spec, "DEV"),
                    L=8,
                    k=np.pi / 4,
                    q=(np.pi / 2, 3 * np.pi / 4),
                    tau=(0.1, 0.27, 0.63, 0.94),
                    g=0.8,
                )
                b = import_r1_binding(rt)
                write_json(output / "binding.json", asdict(b))
                reference = complex(np.asarray(interpret(plan.graph, binder.bind(rt))["F"]).item())
                write_json(
                    output / "demo_result.json",
                    {
                        "case": spec.spec_id,
                        "r1": [reference.real, reference.imag],
                        "compiled": [evaluator(b)["F"].real, evaluator(b)["F"].imag],
                        "ir_file": "r1_scalar_n2.json",
                        "binding_file": "binding.json",
                    },
                )
    write_csv(output / "differential.csv", rows)
    write_csv(output / "graph_stats.csv", stats)
    write_json(
        output / "R1_ORACLE.json",
        {
            "physics_object": "PF_G_FIXED_TIME_ORDER_G2n",
            "graph_object": "SHARED_X_GROUP",
            "atol": 1e-12,
            "rtol": 1e-10,
            "numeric_oracle": "r1_symbolic_build/r1_twoband_build + ReferenceBinder + exec.interpret",
            "independent_oracle": "index_oracle.group_sum_scalar/group_sum_matrix_twoband",
            "frontend": "PlanSpec structure adapter; not reverse engineering of arbitrary Graph objects",
            "primary_lowering": "explicit diagram factor chains with ordered structural hash-consing",
            "recovered_dp": "manual transcription of the R1 recurrence; not automatic CSE rediscovery",
            "validated_orders": {"scalar": [1, 2, 3, 4], "twoband": [1, 2, 3]},
            "rule_limit": "1 <= n <= 6; the admitted orders above the validated ones are not validated",
            "oracle_source_sha256": {
                name: hashlib.sha256((ROOT / name).read_bytes()).hexdigest() for name in ORACLE_SOURCES
            },
            "not_covered": [
                "time integrals",
                "external propagators",
                "the native FEP-DMC weight",
                "variable order",
            ],
        },
    )
    max_error = max(r["max_error_including_Op"] for r in rows)
    return finish(
        "stage 1 (exact compiler)",
        output,
        {
            "ir_json_round_trip": round_trip_ok,
            "all_cases_agree_with_both_oracles": all(r["status"] == "PASS" for r in rows),
            "number_of_cases_matches_the_record": len(rows) == golden()["differential_cases"],
        },
        differential_cases=len(rows),
        max_error_including_Op=max_error,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--output", type=Path, default=ROOT / "build/validation/stage1_exact")
    sys.exit(main(parser.parse_args().output))
