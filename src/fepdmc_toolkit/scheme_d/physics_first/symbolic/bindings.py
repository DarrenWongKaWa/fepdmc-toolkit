"""Deterministic rebind fixtures. Seeds written before any numerical eval."""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np

from .guards import RuntimeBinding, torus_point
from .plan_spec import PREREGISTERED_SCALAR, PlanSpec

NOTES_CASES = Path(__file__).with_name("rebind_cases.csv")
BASE_DEV = 202609120
BASE_RES = 202609130
DEFAULT_L = 4


def _ordered_tau(rng: np.random.Generator, n_tau: int, lo: float = 0.08, hi: float = 1.32) -> tuple[float, ...]:
    raw = np.sort(rng.uniform(lo, hi, size=n_tau))
    for i in range(1, n_tau):
        if raw[i] <= raw[i - 1]:
            raw[i] = raw[i - 1] + 1e-7
    return tuple(float(x) for x in raw)


def base_binding(spec: PlanSpec, split: str) -> RuntimeBinding:
    n = spec.n
    L = DEFAULT_L
    seed = (BASE_DEV if split == "DEV" else BASE_RES) + n
    rng = np.random.default_rng(seed)
    k_index = 0
    q_index = tuple(int(rng.integers(0, L)) for _ in range(n))
    tau = _ordered_tau(rng, 2 * n)
    return RuntimeBinding(
        k=torus_point(k_index, L),
        q=tuple(torus_point(i, L) for i in q_index),
        tau=tau,
        t=1.0,
        omega=0.8,
        g=0.45,
        L=L,
        k_index=k_index,
        q_index=q_index,
    )


def apply_scenario(base: RuntimeBinding, spec: PlanSpec, scenario: str, fixture_seed: int) -> RuntimeBinding:
    rng = np.random.default_rng(int(fixture_seed))
    n = spec.n
    L = int(base.L)
    k_index = int(base.k_index)
    q_index = list(base.q_index)
    tau = list(base.tau)
    t, omega, g = base.t, base.omega, base.g
    if scenario == "change_k":
        k_index = (k_index + 1 + int(rng.integers(0, max(L - 1, 1)))) % L
        if k_index == base.k_index and L > 1:
            k_index = (k_index + 1) % L
    elif scenario == "change_q":
        q_index = [int(rng.integers(0, L)) for _ in range(n)]
        if tuple(q_index) == base.q_index:
            q_index[0] = (q_index[0] + 1) % L
    elif scenario == "change_tau":
        tau = list(_ordered_tau(rng, 2 * n, 0.05, 1.25))
    elif scenario == "change_t":
        t = 1.0 + 0.15 * (1 + int(rng.integers(0, 5)))
    elif scenario == "change_omega":
        omega = 0.5 + 0.1 * (1 + int(rng.integers(0, 6)))
    elif scenario == "change_g":
        g = 0.2 + 0.1 * (1 + int(rng.integers(0, 8)))
        if abs(g - base.g) < 1e-15:
            g = base.g + 0.25
    elif scenario == "g_zero":
        g = 0.0
    elif scenario == "joint_parameters":
        t = 1.25
        omega = 1.05
        g = 0.7
    elif scenario == "joint_sample":
        k_index = (k_index + 2) % L
        q_index = [(qi + 1 + i) % L for i, qi in enumerate(q_index)]
        tau = list(_ordered_tau(rng, 2 * n, 0.12, 1.18))
    elif scenario == "repeat_q_then_separate":
        if n >= 2:
            q_index = [q_index[0]] * n
        else:
            q_index = [(q_index[0] + 1) % L]
    elif scenario == "repeat_dtau_then_separate":
        step = 0.17
        tau = [0.11 + step * i for i in range(2 * n)]
    elif scenario == "change_L_analytic":
        L = 8
        k_index = k_index % L
        q_index = [qi % L for qi in q_index]
    elif scenario == "small_nonzero_g":
        g = 1e-8
    elif scenario == "longer_legal_intervals":
        tau = list(_ordered_tau(rng, 2 * n, 0.2, 2.4))
    elif scenario == "same_input_repeat":
        pass
    elif scenario == "interleaved_contexts":
        k_index = (k_index + 1) % L
        g = 0.55
    else:
        raise ValueError(scenario)
    return RuntimeBinding(
        k=torus_point(k_index, L),
        q=tuple(torus_point(i, L) for i in q_index),
        tau=tuple(float(x) for x in tau),
        t=float(t),
        omega=float(omega),
        g=float(g),
        L=int(L),
        k_index=int(k_index),
        q_index=tuple(int(i) for i in q_index),
        material_id=base.material_id,
        gauge_id=base.gauge_id,
        data_version=base.data_version,
        finite_model=True,
        delta=base.delta,
        gap=base.gap,
    )


def load_rebind_cases(path: Path = NOTES_CASES) -> list[dict]:
    return list(csv.DictReader(path.open()))


def materialize_case(row: dict) -> tuple[PlanSpec, RuntimeBinding]:
    spec = PREREGISTERED_SCALAR[row["spec_id"]]
    split = row["split"]
    base = base_binding(spec, split)
    binding = apply_scenario(base, spec, row["scenario"], int(row["fixture_seed"]))
    return spec, binding


def write_seeded_bindings(out_csv: Path, path: Path = NOTES_CASES) -> list[dict]:
    rows = []
    for row in load_rebind_cases(path):
        spec, binding = materialize_case(row)
        rec = {
            "case_id": row["case_id"],
            "spec_id": spec.spec_id,
            "split": row["split"],
            "scenario": row["scenario"],
            "fixture_seed": row["fixture_seed"],
            "n": spec.n,
            "k": repr(binding.k),
            "q": ";".join(repr(x) for x in binding.q),
            "tau": ";".join(repr(x) for x in binding.tau),
            "t": repr(binding.t),
            "omega": repr(binding.omega),
            "g": repr(binding.g),
            "L": binding.L,
            "finite_model": binding.finite_model,
            "state": "SEEDED_NO_EVAL",
        }
        rows.append(rec)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with out_csv.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    return rows
