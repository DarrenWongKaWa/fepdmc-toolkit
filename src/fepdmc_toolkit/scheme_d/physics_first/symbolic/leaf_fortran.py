"""Unique writer: scalar leaf table -> Fortran assignments from raw physics inputs.

R1 and CSE must call this, not a second formula copy. Coefficients 2 and 1 in
xi=2 t (1-cos p) are the frozen dispersion, not sample values.
"""

from __future__ import annotations

from .leaves import ElectronProp, LeafTable, PhononProp, PrefactorProp, VertexProp
from .plan_spec import PlanSpec


def emit_leaf_assignments(
    leaves: LeafTable,
    spec: PlanSpec,
    variables: dict,
    *,
    include_l_guard: bool = True,
) -> list[str]:
    if spec.electron_interface != "scalar" or spec.bands != 1:
        raise ValueError("fused leaf Fortran is scalar Holstein only")
    lines: list[str] = []
    for name, item in variables.items():
        desc = leaves.descriptors.get(name)
        if desc is None:
            raise KeyError(f"DAG input {name} is not in the leaf table")
        off = int(item["offset"]) + 1
        size = int(item["size"])
        if size != 1:
            raise ValueError(f"scalar fused kernel expects size-1 leaf {name}")
        lines.extend(_assign(desc, off, include_l_guard=include_l_guard))
    return lines


def _assign(desc, off: int, *, include_l_guard: bool = True) -> list[str]:
    dest = f"x({off})"
    if isinstance(desc, PrefactorProp):
        if desc.normalization_spec != "g_over_sqrt_L":
            raise ValueError(desc.normalization_spec)
        p = 2 * int(desc.n)
        guard = (
            [
                "if (Lval <= 0.0d0) then",
                "status=3",
                "return",
                "end if",
            ]
            if include_l_guard
            else []
        )
        return guard + [
            "gv = g / sqrt(Lval)",
            f"{dest} = cmplx(gv ** {p}, 0.0d0, kind=c_double)",
        ]
    if isinstance(desc, ElectronProp):
        if desc.electron_interface != "scalar":
            raise ValueError(desc.electron_interface)
        lines = [f"k_eff = {desc.momentum.k_coeff}.0d0 * k"]
        for slot, coeff in desc.momentum.q_coeffs:
            lines.append(f"k_eff = k_eff + ({coeff}.0d0) * q({int(slot) + 1})")
        lines.append(f"dtau = tau({desc.tau_right_idx + 1}) - tau({desc.tau_left_idx + 1})")
        lines.append("xi = 2.0d0 * t * (1.0d0 - cos(k_eff))")
        lines.append(f"{dest} = cmplx(exp(-xi * dtau), 0.0d0, kind=c_double)")
        return lines
    if isinstance(desc, PhononProp):
        if desc.layout_bands != 1:
            raise ValueError("scalar fused kernel")
        return [
            f"dtau = tau({desc.tau_close_idx + 1}) - tau({desc.tau_open_idx + 1})",
            f"{dest} = cmplx(exp(-omega * dtau), 0.0d0, kind=c_double)",
        ]
    if isinstance(desc, VertexProp):
        raise ValueError("scalar fused kernel has no separate vertex leaf")
    raise TypeError(type(desc))


def leaf_fortran_text(leaves: LeafTable, spec: PlanSpec, variables: dict) -> str:
    return "\n".join(emit_leaf_assignments(leaves, spec, variables)) + "\n"
