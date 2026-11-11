"""Independent chord enumeration, canonical factors, and bounded IR validation."""

from __future__ import annotations

import json

from .ir import (
    KIND_ELECTRON,
    KIND_PHONON,
    KIND_PREFACTOR,
    KIND_VERTEX,
    MODEL_SCALAR,
    MODEL_TWOBAND,
    SERIES_FULL,
    Diagram,
    DiagramIR,
    Factor,
    Family,
    MomentumForm,
    Variables,
)

Pairing = tuple[tuple[int, int], ...]


class IrError(ValueError):
    """Unsupported or inconsistent diagram semantics."""


def enumerate_pairings(n: int) -> tuple[Pairing, ...]:
    """Enumerate perfect matchings by independent leftmost-vertex recursion."""
    if type(n) is not int or not 1 <= n <= 6:
        raise ValueError("supported resource bound: integer 1 <= n <= 6")

    def rec(pts: tuple[int, ...]) -> list[Pairing]:
        if not pts:
            return [()]
        out: list[Pairing] = []
        for i in range(1, len(pts)):
            for tail in rec(pts[1:i] + pts[i + 1 :]):
                out.append(tuple(sorted(((pts[0], pts[i]),) + tail)))
        return out

    return tuple(rec(tuple(range(1, 2 * n + 1))))


def chord_slots(pairing: Pairing) -> dict[tuple[int, int], int]:
    return {edge: i for i, edge in enumerate(sorted(pairing))}


def open_chord_slots(pairing: Pairing, cut: int) -> tuple[int, ...]:
    slots = chord_slots(pairing)
    return tuple(sorted(slots[e] for e in pairing if e[0] <= cut < e[1]))


def diagram_id(pairing: Pairing) -> str:
    return "C_" + "_".join(f"({a},{b})" for a, b in sorted(pairing))


def _electron_factor(pairing: Pairing, v: int) -> Factor:
    slots = open_chord_slots(pairing, v)
    return Factor(
        KIND_ELECTRON,
        interval=(v - 1, v),
        open_chord_slots=slots,
        k_form=MomentumForm(1, tuple((s, -1) for s in slots)),
    )


def build_scalar_diagram(pairing: Pairing) -> Diagram:
    n = len(pairing)
    factors = [_electron_factor(pairing, v) for v in range(1, 2 * n)]
    for edge, slot in chord_slots(pairing).items():
        factors.append(
            Factor(
                KIND_PHONON,
                chord_slot=slot,
                tau_open_idx=edge[0] - 1,
                tau_close_idx=edge[1] - 1,
            )
        )
    return Diagram(diagram_id(pairing), tuple(sorted(pairing)), 1, 1, tuple(factors))


def build_twoband_diagram(pairing: Pairing) -> Diagram:
    n = len(pairing)
    slots = chord_slots(pairing)
    factors: list[Factor] = []
    for v in range(1, 2 * n + 1):
        for edge in sorted(pairing):
            if v in edge:
                factors.append(
                    Factor(
                        KIND_VERTEX,
                        vertex=v,
                        chord_slot=slots[edge],
                        direction="emit" if v == edge[0] else "absorb",
                    )
                )
        if v < 2 * n:
            factors.append(_electron_factor(pairing, v))
    for edge, slot in slots.items():
        factors.append(
            Factor(
                KIND_PHONON,
                chord_slot=slot,
                tau_open_idx=edge[0] - 1,
                tau_close_idx=edge[1] - 1,
            )
        )
    return Diagram(diagram_id(pairing), tuple(sorted(pairing)), 1, 1, tuple(factors))


def scalar_family(n: int) -> Family:
    return Family(
        family_id=f"holstein_scalar_vacuum_fixed_time_n{n}",
        model=MODEL_SCALAR,
        order=n,
        bands=1,
        electron_interface="scalar",
        propagator_convention="physical_unshifted",
        normalization_spec="g_over_sqrt_L",
        series=SERIES_FULL,
        dispersion="xi_k=2t(1-cos k)",
        phonon_spec="einstein_vacuum_T0",
        rule_version="diagram-compiler-v1",
        global_factors=(Factor(KIND_PREFACTOR, normalization="g_over_sqrt_L", exponent="2n"),),
    )


def twoband_family(n: int) -> Family:
    return Family(
        family_id=f"holstein_twoband_hermitian_fixed_time_n{n}",
        model=MODEL_TWOBAND,
        order=n,
        bands=2,
        electron_interface="matrix2",
        propagator_convention="physical_unshifted",
        normalization_spec="vertex_contains_g_over_sqrt_L",
        series=SERIES_FULL,
        dispersion="xi_k=2t(1-cos k)",
        phonon_spec="einstein_vacuum_T0",
        rule_version="diagram-compiler-v1",
        global_factors=(),
    )


def _variables(family: Family) -> Variables:
    params = ("t", "omega", "g", "L")
    if family.model == MODEL_TWOBAND:
        params += ("delta", "gap")
    return Variables(
        external=("k",),
        chord_slots=tuple(f"q[{i}]" for i in range(family.order)),
        chord_slot_rule="first_opening_rank",
        time_slots=tuple(f"tau[{i}]" for i in range(2 * family.order)),
        time_rule="vertex_times_strictly_increasing",
        params=params,
    )


def build_scalar_ir(n: int) -> DiagramIR:
    pairings = enumerate_pairings(n)
    family = scalar_family(n)
    return DiagramIR(family, _variables(family), tuple(build_scalar_diagram(p) for p in pairings))


def build_twoband_ir(n: int) -> DiagramIR:
    pairings = enumerate_pairings(n)
    family = twoband_family(n)
    return DiagramIR(family, _variables(family), tuple(build_twoband_diagram(p) for p in pairings))


def _require(cond: bool, msg: str) -> None:
    if not cond:
        raise IrError(msg)


def validate_ir(ir: DiagramIR) -> None:
    family = ir.family
    n = family.order
    _require(ir.object_type == "SHARED_X_GROUP", "only SHARED_X_GROUP is supported")
    _require(ir.schema_version == "1.0.0", f"unsupported schema {ir.schema_version}")
    _require(type(n) is int and 1 <= n <= 6, "supported resource bound: integer 1 <= n <= 6")
    _require(type(family.bands) is int, "bands must be an integer")
    _require(family.model in (MODEL_SCALAR, MODEL_TWOBAND), f"unsupported model {family.model}")
    _require(isinstance(family.family_id, str) and bool(family.family_id), "family_id must be a label")
    expected_family = scalar_family(n) if family.model == MODEL_SCALAR else twoband_family(n)
    for field in Family.__slots__:
        if field != "family_id":
            _require(getattr(family, field) == getattr(expected_family, field), f"unsupported family {field}")
    _require(ir.variables == _variables(expected_family), "variable binding contract mismatch")
    pairings = enumerate_pairings(n)
    _require(tuple(d.pairing for d in ir.diagrams) == pairings, "FULL series requires canonical pairing order")
    builder = build_scalar_diagram if family.model == MODEL_SCALAR else build_twoband_diagram
    for diagram, pairing in zip(ir.diagrams, pairings, strict=True):
        actual = json.dumps(diagram.to_dict(), sort_keys=True, allow_nan=False)
        expected = json.dumps(builder(pairing).to_dict(), sort_keys=True)
        _require(actual == expected, f"diagram factors or identity mismatch: {diagram.diagram_id}")
