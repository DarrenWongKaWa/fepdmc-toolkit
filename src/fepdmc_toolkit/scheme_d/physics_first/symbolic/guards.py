"""Runtime binding guards. Reject before kernel. Do not silently sort or snap."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .plan_spec import PlanSpec


class GuardError(ValueError):
    def __init__(self, code: str, message: str) -> None:
        self.code = code
        super().__init__(f"{code}: {message}")


@dataclass(frozen=True)
class RuntimeBinding:
    k: float
    q: tuple[float, ...]
    tau: tuple[float, ...]
    t: float
    omega: float
    g: float
    L: int
    delta: float = 0.35
    gap: float = 0.5
    material_id: str = "synthetic"
    gauge_id: str = "default"
    data_version: str = "analytic_holstein_v1"
    finite_model: bool = True
    continuous_extension: bool = False
    k_index: int = 0
    q_index: tuple[int, ...] = ()

    def q_array(self) -> np.ndarray:
        return np.asarray(self.q, dtype=np.float64)

    def tau_array(self) -> np.ndarray:
        return np.asarray(self.tau, dtype=np.float64)


def torus_point(index: int, L: int) -> float:
    return float(2.0 * np.pi * (int(index) % L) / L)


def is_on_torus(value: float, L: int, atol: float = 1e-12) -> bool:
    if L <= 0:
        return False
    idx = int(np.round(value * L / (2.0 * np.pi)))
    snapped = 2.0 * np.pi * (idx % L) / L
    wrapped = 2.0 * np.pi * ((idx % L) + L) / L
    return bool(min(abs(value - snapped), abs(value - snapped - 2.0 * np.pi), abs(value - wrapped)) <= atol)


def _finite(*vals: float) -> bool:
    return all(np.isfinite(v) for v in vals)


def validate_binding(spec: PlanSpec, binding: RuntimeBinding) -> None:
    q = binding.q_array()
    tau = binding.tau_array()
    if not isinstance(binding.L, (int, np.integer)) or int(binding.L) <= 0:
        raise GuardError("G04", f"invalid_positive_L L={binding.L}")
    if q.shape != (spec.n,):
        raise GuardError("G02", f"wrong_q_slot_count got {q.shape} want ({spec.n},)")
    if tau.shape != (2 * spec.n,):
        raise GuardError("G02", f"wrong_tau_count got {tau.shape} want ({2 * spec.n},)")
    fields = [binding.k, binding.t, binding.omega, binding.g, float(binding.L), binding.delta, binding.gap]
    fields.extend(float(x) for x in q)
    fields.extend(float(x) for x in tau)
    if not _finite(*fields):
        raise GuardError("G03", "nonfinite_input")
    if spec.bands == 1 and spec.electron_interface != "scalar":
        raise GuardError("G07", "scalar spec requires scalar interface")
    if spec.bands == 2 and spec.electron_interface != "matrix2":
        raise GuardError("G07", "two-band spec requires matrix2 interface")
    if np.any(np.diff(tau) <= 0):
        raise GuardError("G01", "unsorted_times; do not silently sort")
    if binding.material_id != spec.expected_material_id or binding.gauge_id != spec.expected_gauge_id:
        raise GuardError("G08", "material_gauge_mismatch")
    if binding.data_version != spec.expected_data_version:
        raise GuardError("G08", "data_version_mismatch")
    if binding.finite_model:
        if binding.continuous_extension:
            raise GuardError("G05", "finite_model cannot also be continuous_extension")
        if not is_on_torus(binding.k, int(binding.L)):
            raise GuardError("G05", "offgrid_k_as_finite_model")
        for i, qi in enumerate(q):
            if not is_on_torus(float(qi), int(binding.L)):
                raise GuardError("G05", f"offgrid_q_as_finite_model slot={i}")
    elif not binding.continuous_extension:
        raise GuardError("G05", "offgrid requires explicit continuous_extension; not a finite-model/Fock point")
    if spec.propagator_convention != "physical_unshifted":
        raise GuardError("G06", "unsupported propagator convention")


def refuse_unit_g_raw(spec: PlanSpec) -> None:
    """B1 native unit-G skip is illegal on this raw physical kernel."""
    from fepdmc_toolkit.scheme_d.rewrites import RewriteRejected

    if spec.propagator_convention == "physical_unshifted":
        raise RewriteRejected("false_unit_G_in_raw_scalar")
