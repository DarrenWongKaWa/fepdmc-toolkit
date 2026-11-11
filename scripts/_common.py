"""Helpers shared by the diagram-compiler validation stages (scripts/dc_stage*.py).

Every stage writes its records to an output directory, ends with ``finish`` and
exits with status 1 if one of its checks fails.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np

from fepdmc_toolkit.diagram_compiler import Binding
from fepdmc_toolkit.diagram_compiler.evaluator import torus_point

ROOT = Path(__file__).resolve().parents[1]
RECORDS = ROOT / "validation" / "diagram_compiler"

# Three-state chain of the detailed-balance checks: target, cheap weight, proposals.
PI = np.array([0.5, 0.3, 0.2])
CHEAP = np.array([0.1, 0.55, 0.35])
Q_SYMMETRIC = np.array([[0.0, 0.5, 0.5], [0.5, 0.0, 0.5], [0.5, 0.5, 0.0]])
Q_ASYMMETRIC = np.array([[0.0, 0.8, 0.2], [0.3, 0.0, 0.7], [0.6, 0.4, 0.0]])
BALANCE_TOL = 1e-12


def write_csv(path: Path, rows: list[dict]) -> None:
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0].keys()), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def write_json(path: Path, payload) -> None:
    path.write_text(json.dumps(payload, indent=2, sort_keys=True, allow_nan=False) + "\n")


def golden() -> dict:
    """Reference values of the natively validated two-band n = 1 kernel."""
    return json.loads((RECORDS / "golden.json").read_text())


def twoband_binding(variant: str) -> Binding:
    """The two sample points x and y of the two-band n = 1 demonstrations."""
    L = 4
    if variant == "x":
        k_idx, q0, tau_step, t, omega, g = 0, 1, 0.05, 1.0, 0.8, 0.5
    else:
        k_idx, q0, tau_step, t, omega, g = 1, 2, 0.07, 1.1, 0.9, 0.4
    return Binding(
        k=torus_point(k_idx, L),
        q=(torus_point(q0, L),),
        tau=(tau_step, 2 * tau_step),
        t=t,
        omega=omega,
        g=g,
        L=L,
        delta=0.35,
        gap=0.5,
    )


def balance_residuals(p: np.ndarray) -> tuple[float, float]:
    """Largest detailed-balance and stationarity residuals of the transition matrix ``p`` for PI."""
    n = len(PI)
    balance = max(abs(PI[i] * p[i, j] - PI[j] * p[j, i]) for i in range(n) for j in range(n))
    return float(balance), float(np.max(np.abs(PI @ p - PI)))


def finish(stage: str, output: Path, checks: dict[str, bool], **numbers) -> int:
    """Write summary.json, print the outcome and return the exit status of the stage."""
    checks = {name: bool(ok) for name, ok in checks.items()}
    failed = [name for name, ok in checks.items() if not ok]
    status = "FAIL" if failed else "PASS"
    write_json(output / "summary.json", {"stage": stage, "status": status, "checks": checks, **numbers})
    print(f"{stage}: {status} ({len(checks) - len(failed)} of {len(checks)} checks), records in {output}")
    for name in failed:
        print(f"  failed: {name}")
    return 1 if failed else 0
