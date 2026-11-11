"""Compiler-oracle fixtures for stage-6 native validation. Not DA mathematics."""

from __future__ import annotations

import math
from dataclasses import dataclass

from .evaluator import Binding, torus_point
from .native_interpret import interpret_native_kernel
from .native_ir import NativeKernelIR


def demo_binding(n: int = 1, *, variant: str = "x") -> Binding:
    L = 4
    if variant == "x":
        k_idx, q0, tau_step, t, omega, g = 0, 1, 0.05, 1.0, 0.8, 0.5
    else:
        k_idx, q0, tau_step, t, omega, g = 1, 2, 0.07, 1.1, 0.9, 0.4
    return Binding(
        k=torus_point(k_idx, L),
        q=tuple(torus_point(q0 + i, L) for i in range(n)),
        tau=tuple(tau_step * (i + 1) for i in range(2 * n)),
        t=t,
        omega=omega,
        g=g,
        L=L,
        delta=0.35,
        gap=0.5,
    )


def binding_to_dict(binding: Binding, layout) -> dict:
    payload = {}
    for field in layout.fields:
        value = getattr(binding, field.name)
        if field.shape == ():
            payload[field.name] = int(value) if field.name == "L" else float(value)
        else:
            payload[field.name] = [float(v) for v in value]
    return payload


def fixture_text(layout, x: Binding, y: Binding, log_q, u1, u2, exact_x_valid, exact_log_weight_x) -> str:
    lines: list[str] = []
    for binding in (x, y):
        for field in layout.fields:
            value = getattr(binding, field.name)
            if field.shape == ():
                lines.append(str(value))
            else:
                lines.append(" ".join(str(float(v)) for v in value))
    lines.append(f"{float(log_q)} {float(u1)} {float(u2)}")
    lines.append(f"{1 if exact_x_valid else 0} {float(exact_log_weight_x)}")
    return "\n".join(lines) + "\n"


@dataclass(frozen=True)
class OracleCase:
    name: str
    text: str
    expected: dict
    x: Binding
    y: Binding
    log_q: float
    u1: float
    u2: float
    exact_x_valid: bool
    exact_log_weight_x: float


def _expected(result) -> dict:
    return {
        "accepted": result.accepted,
        "ell_R": result.ell_R,
        "ell_R_valid": result.ell_R is not None,
        "ell_hat": result.ell_hat,
        "exact_y_evaluated": result.exact_y_evaluated,
        "reject_stage": result.reject_stage,
        "status": result.status,
    }


def _case(nir: NativeKernelIR, name, x, y, log_q, u1, u2, exact_x_valid=False, exact_log_weight_x=0.0) -> OracleCase:
    py = interpret_native_kernel(
        nir,
        x=x,
        y=y,
        log_q=log_q,
        u1=u1,
        u2=u2,
        exact_x_valid=exact_x_valid,
        exact_log_weight_x=exact_log_weight_x,
    )
    return OracleCase(
        name=name,
        text=fixture_text(nir.layout, x, y, log_q, u1, u2, exact_x_valid, exact_log_weight_x),
        expected=_expected(py),
        x=x,
        y=y,
        log_q=log_q,
        u1=u1,
        u2=u2,
        exact_x_valid=exact_x_valid,
        exact_log_weight_x=exact_log_weight_x,
    )


def twoband_oracle_cases(nir: NativeKernelIR, kernel) -> list[OracleCase]:
    x, y = demo_binding(1, variant="x"), demo_binding(1, variant="y")
    probe = interpret_native_kernel(nir, x=x, y=y, log_q=0.0, u1=1e-16, u2=1e-16)
    alpha2 = math.exp(min(0.0, probe.ell_R - probe.ell_hat))
    u2_reject = (alpha2 + 1.0) / 2.0
    log_px = kernel.exact_log_weight(x)
    zero = Binding(k=x.k, q=x.q, tau=x.tau, t=x.t, omega=x.omega, g=0.0, L=x.L, delta=x.delta, gap=x.gap)
    return [
        _case(nir, "stage1_reject", x, y, 0.0, 0.9, 0.1),
        _case(nir, "stage2_reject", x, y, 0.0, 1e-16, u2_reject),
        _case(nir, "accept", x, y, 0.0, 1e-16, 1e-16),
        _case(nir, "asymmetric_q", x, y, 0.4, 1e-16, 1e-16),
        _case(nir, "cached_exact_x", x, y, 0.0, 1e-16, 1e-16, True, log_px),
        _case(nir, "uncached_exact_x", x, y, 0.0, 1e-16, 1e-16, False, 0.0),
        _case(nir, "invalid_target_g0", zero, y, 0.0, 0.5, 0.5),
        _case(nir, "invalid_proposal_ratio", x, y, float("inf"), 0.5, 0.5),
        _case(nir, "invalid_u1", x, y, 0.0, -0.1, 0.5),
        _case(nir, "invalid_u2", x, y, 0.0, 0.5, 1.1),
    ]
