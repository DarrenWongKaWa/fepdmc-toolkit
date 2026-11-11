"""Numeric identity gate. Near-zero uses absolute, not pure relative error."""

from __future__ import annotations

import math

ABS_ATOL = 1e-12
REL_RTOL = 1e-10


def reference_scale(reference: complex | float) -> float:
    return abs(complex(reference))


def numeric_gate(value: complex | float, reference: complex | float) -> bool:
    if not all(math.isfinite(x) for z in (complex(value), complex(reference)) for x in (z.real, z.imag)):
        return False
    err = abs(complex(value) - complex(reference))
    scale = reference_scale(reference)
    # A fixed 1e-12 floor cannot resolve an amplitude of 1e-223. Require its
    # relative value, or compare normalized matrices/log-scales instead.
    tol = REL_RTOL * scale if 0 < scale < ABS_ATOL else ABS_ATOL + REL_RTOL * scale
    return err <= tol


def gate_margin(value: complex | float, reference: complex | float) -> float:
    return abs(complex(value) - complex(reference)) - (ABS_ATOL + REL_RTOL * reference_scale(reference))


def require_finite(z: complex | float, name: str = "value") -> complex:
    c = complex(z)
    if not math.isfinite(c.real) or not math.isfinite(c.imag):
        raise ValueError(f"non-finite {name}")
    return c
