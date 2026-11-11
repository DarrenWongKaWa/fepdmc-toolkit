"""Explicit DiagramIR for the R1 fixed-order grouped-diagram family.

This package is the explicit diagram-level representation that the R1
evaluator (fepdmc_toolkit.scheme_d, a compile-once recursive implementation)
uses only implicitly.

Three layers:

A. PHYSICAL SEMANTICS  -- ir.py (this module): family, variables, diagrams,
   per-diagram semantic factors. Understandable and serializable alone.
B. EVALUATION SEMANTICS -- graph.py + evaluator.py: lowering rules and
   the exact evaluation DAG + primitives.
C. OPTIMIZATION METADATA -- canonical.py structural keys; sharing decisions
   live in the lowering, never in the IR itself.

Scope: exactly the two families R1 supports,
Holstein scalar vacuum n>=1 and Holstein two-band hermitian n>=1, fixed
ordered vertex times, series FULL_ALL_PAIRINGS. Not an arbitrary-diagram
compiler and not a FEP-DMC reproduction.
"""

from __future__ import annotations

from dataclasses import dataclass
from numbers import Integral

SCHEMA_VERSION = "1.0.0"

OBJECT_SHARED_X_GROUP = "SHARED_X_GROUP"
RULE_VERSION = "diagram-compiler-v1"

MODEL_SCALAR = "Holstein_scalar_vacuum"
MODEL_TWOBAND = "Holstein_twoband_hermitian"
SERIES_FULL = "FULL_ALL_PAIRINGS"

KIND_PREFACTOR = "prefactor"
KIND_ELECTRON = "electron_propagator"
KIND_PHONON = "phonon_propagator"
KIND_VERTEX = "vertex"


def _rev(t):
    if isinstance(t, tuple):
        return [_rev(x) for x in t]
    if isinstance(t, MomentumForm):
        return t.to_dict()
    if isinstance(t, Factor):
        return t.to_dict()
    return t


def _integer(value) -> int:
    if isinstance(value, bool) or not isinstance(value, Integral):
        raise ValueError(f"structural index must be an integer, got {value!r}")
    return int(value)


@dataclass(frozen=True, slots=True, init=False)
class MomentumForm:
    """Immutable integer form k_coeff*k + sum(coeff*q[slot])."""

    k_coeff: int
    q_terms: tuple[tuple[int, int], ...]

    def __init__(self, k_coeff: int, q_terms: tuple[tuple[int, int], ...]) -> None:
        object.__setattr__(self, "k_coeff", _integer(k_coeff))
        object.__setattr__(self, "q_terms", tuple((_integer(s), _integer(c)) for s, c in q_terms))

    def to_dict(self) -> dict:
        return {"k_coeff": self.k_coeff, "q_terms": [[s, c] for s, c in self.q_terms]}

    @classmethod
    def from_dict(cls, d: dict) -> MomentumForm:
        return cls(d["k_coeff"], tuple((s, c) for s, c in d["q_terms"]))

    def _tup(self) -> tuple:
        return (self.k_coeff, self.q_terms)

    def __eq__(self, other: object) -> bool:
        return isinstance(other, MomentumForm) and self._tup() == other._tup()

    def __repr__(self) -> str:
        if not self.q_terms:
            return f"k*{self.k_coeff}"
        return f"k*{self.k_coeff} " + " ".join(f"q[{s}]*{c}" for s, c in self.q_terms)


class Factor:
    """One semantic diagram factor. Fields unused by `kind` stay None.

    kinds:
      prefactor          normalization, exponent
      electron_propagator interval, open_chord_slots, k_form
      phonon_propagator   chord_slot, tau_open_idx, tau_close_idx
      vertex              vertex (1-based), chord_slot, direction
    """

    __slots__ = (
        "kind",
        "normalization",
        "exponent",
        "interval",
        "open_chord_slots",
        "k_form",
        "chord_slot",
        "tau_open_idx",
        "tau_close_idx",
        "vertex",
        "direction",
    )

    def __init__(
        self,
        kind: str,
        *,
        normalization: str | None = None,
        exponent: str | None = None,
        interval: tuple[int, int] | None = None,
        open_chord_slots: tuple[int, ...] | None = None,
        k_form: MomentumForm | None = None,
        chord_slot: int | None = None,
        tau_open_idx: int | None = None,
        tau_close_idx: int | None = None,
        vertex: int | None = None,
        direction: str | None = None,
    ) -> None:
        self.kind = kind
        self.normalization = normalization
        self.exponent = exponent
        self.interval = None if interval is None else (_integer(interval[0]), _integer(interval[1]))
        self.open_chord_slots = None if open_chord_slots is None else tuple(_integer(s) for s in open_chord_slots)
        self.k_form = k_form
        self.chord_slot = None if chord_slot is None else _integer(chord_slot)
        self.tau_open_idx = None if tau_open_idx is None else _integer(tau_open_idx)
        self.tau_close_idx = None if tau_close_idx is None else _integer(tau_close_idx)
        self.vertex = None if vertex is None else _integer(vertex)
        self.direction = direction

    def _tup(self) -> tuple:
        return (
            self.kind,
            self.normalization,
            self.exponent,
            self.interval,
            self.open_chord_slots,
            self.k_form,
            self.chord_slot,
            self.tau_open_idx,
            self.tau_close_idx,
            self.vertex,
            self.direction,
        )

    def __eq__(self, other: object) -> bool:
        return isinstance(other, Factor) and self._tup() == other._tup()

    def __repr__(self) -> str:
        fields = {s: getattr(self, s) for s in self.__slots__ if getattr(self, s) is not None}
        inner = ",".join(f"{k}={v}" for k, v in fields.items() if k != "kind")
        return f"{self.kind}({inner})"

    def to_dict(self) -> dict:
        d: dict = {"kind": self.kind}
        if self.normalization is not None:
            d["normalization"] = self.normalization
        if self.exponent is not None:
            d["exponent"] = self.exponent
        if self.interval is not None:
            d["interval"] = list(self.interval)
        if self.open_chord_slots is not None:
            d["open_chord_slots"] = list(self.open_chord_slots)
        if self.k_form is not None:
            d["k_form"] = self.k_form.to_dict()
        if self.chord_slot is not None:
            d["chord_slot"] = self.chord_slot
        if self.tau_open_idx is not None:
            d["tau_open_idx"] = self.tau_open_idx
        if self.tau_close_idx is not None:
            d["tau_close_idx"] = self.tau_close_idx
        if self.vertex is not None:
            d["vertex"] = self.vertex
        if self.direction is not None:
            d["direction"] = self.direction
        return d

    @classmethod
    def from_dict(cls, d: dict) -> Factor:
        return cls(
            d["kind"],
            normalization=d.get("normalization"),
            exponent=d.get("exponent"),
            interval=tuple(d["interval"]) if "interval" in d else None,
            open_chord_slots=tuple(d["open_chord_slots"]) if "open_chord_slots" in d else None,
            k_form=MomentumForm.from_dict(d["k_form"]) if "k_form" in d else None,
            chord_slot=d.get("chord_slot"),
            tau_open_idx=d.get("tau_open_idx"),
            tau_close_idx=d.get("tau_close_idx"),
            vertex=d.get("vertex"),
            direction=d.get("direction"),
        )


class Diagram:
    """One diagram C: chord pairing identity plus explicit semantic factors."""

    __slots__ = ("diagram_id", "pairing", "sign", "multiplicity", "factors")

    def __init__(
        self,
        diagram_id: str,
        pairing: tuple[tuple[int, int], ...],
        sign: int,
        multiplicity: int,
        factors: tuple[Factor, ...],
    ) -> None:
        self.diagram_id = str(diagram_id)
        self.pairing = tuple((_integer(a), _integer(b)) for a, b in pairing)
        self.sign = _integer(sign)
        self.multiplicity = _integer(multiplicity)
        self.factors = tuple(factors)

    def __eq__(self, other: object) -> bool:
        return isinstance(other, Diagram) and self._tup() == other._tup()

    def _tup(self) -> tuple:
        return (self.diagram_id, self.pairing, self.sign, self.multiplicity, self.factors)

    def __repr__(self) -> str:
        return f"{self.diagram_id} pairing={self.pairing} sign={self.sign} mult={self.multiplicity}"

    def to_dict(self) -> dict:
        return {
            "diagram_id": self.diagram_id,
            "pairing": [list(e) for e in self.pairing],
            "sign": self.sign,
            "multiplicity": self.multiplicity,
            "factors": [f.to_dict() for f in self.factors],
        }

    @classmethod
    def from_dict(cls, d: dict) -> Diagram:
        return cls(
            d["diagram_id"],
            tuple((a, b) for a, b in d["pairing"]),
            d["sign"],
            d["multiplicity"],
            tuple(Factor.from_dict(f) for f in d["factors"]),
        )


class Family:
    """The supported diagram family: model conventions + global factors."""

    __slots__ = (
        "family_id",
        "model",
        "order",
        "bands",
        "electron_interface",
        "propagator_convention",
        "normalization_spec",
        "series",
        "dispersion",
        "phonon_spec",
        "rule_version",
        "global_factors",
    )

    def __init__(
        self,
        family_id: str,
        model: str,
        order: int,
        bands: int,
        electron_interface: str,
        propagator_convention: str,
        normalization_spec: str,
        series: str,
        dispersion: str,
        phonon_spec: str,
        rule_version: str,
        global_factors: tuple[Factor, ...],
    ) -> None:
        self.family_id = family_id
        self.model = model
        self.order = _integer(order)
        self.bands = _integer(bands)
        self.electron_interface = electron_interface
        self.propagator_convention = propagator_convention
        self.normalization_spec = normalization_spec
        self.series = series
        self.dispersion = dispersion
        self.phonon_spec = phonon_spec
        self.rule_version = rule_version
        self.global_factors = tuple(global_factors)

    def __eq__(self, other: object) -> bool:
        return isinstance(other, Family) and self._tup() == other._tup()

    def _tup(self) -> tuple:
        return (
            self.family_id,
            self.model,
            self.order,
            self.bands,
            self.electron_interface,
            self.propagator_convention,
            self.normalization_spec,
            self.series,
            self.dispersion,
            self.phonon_spec,
            self.rule_version,
            self.global_factors,
        )

    def to_dict(self) -> dict:
        return {
            "family_id": self.family_id,
            "model": self.model,
            "order": self.order,
            "bands": self.bands,
            "electron_interface": self.electron_interface,
            "propagator_convention": self.propagator_convention,
            "normalization_spec": self.normalization_spec,
            "series": self.series,
            "dispersion": self.dispersion,
            "phonon_spec": self.phonon_spec,
            "rule_version": self.rule_version,
            "global_factors": [f.to_dict() for f in self.global_factors],
        }

    @classmethod
    def from_dict(cls, d: dict) -> Family:
        return cls(
            d["family_id"],
            d["model"],
            d["order"],
            d["bands"],
            d["electron_interface"],
            d["propagator_convention"],
            d["normalization_spec"],
            d["series"],
            d["dispersion"],
            d["phonon_spec"],
            d["rule_version"],
            tuple(Factor.from_dict(f) for f in d["global_factors"]),
        )


class Variables:
    """Variable scopes with the binding rules made explicit."""

    __slots__ = (
        "external",
        "chord_slots",
        "chord_slot_rule",
        "time_slots",
        "time_rule",
        "params",
    )

    def __init__(
        self,
        external: tuple[str, ...],
        chord_slots: tuple[str, ...],
        chord_slot_rule: str,
        time_slots: tuple[str, ...],
        time_rule: str,
        params: tuple[str, ...],
    ) -> None:
        self.external = tuple(external)
        self.chord_slots = tuple(chord_slots)
        self.chord_slot_rule = chord_slot_rule
        self.time_slots = tuple(time_slots)
        self.time_rule = time_rule
        self.params = tuple(params)

    def __eq__(self, other: object) -> bool:
        return isinstance(other, Variables) and self._tup() == other._tup()

    def _tup(self) -> tuple:
        return (
            self.external,
            self.chord_slots,
            self.chord_slot_rule,
            self.time_slots,
            self.time_rule,
            self.params,
        )

    def to_dict(self) -> dict:
        return {
            "external": list(self.external),
            "chord_slots": list(self.chord_slots),
            "chord_slot_rule": self.chord_slot_rule,
            "time_slots": list(self.time_slots),
            "time_rule": self.time_rule,
            "params": list(self.params),
        }

    @classmethod
    def from_dict(cls, d: dict) -> Variables:
        return cls(
            tuple(d["external"]),
            tuple(d["chord_slots"]),
            d["chord_slot_rule"],
            tuple(d["time_slots"]),
            d["time_rule"],
            tuple(d["params"]),
        )


class DiagramIR:
    """The full explicit IR for one supported family instance."""

    __slots__ = ("schema_version", "object_type", "family", "variables", "diagrams")

    def __init__(
        self,
        family: Family,
        variables: Variables,
        diagrams: tuple[Diagram, ...],
        schema_version: str = SCHEMA_VERSION,
        object_type: str = OBJECT_SHARED_X_GROUP,
    ) -> None:
        self.schema_version = schema_version
        self.object_type = object_type
        self.family = family
        self.variables = variables
        self.diagrams = tuple(diagrams)

    def __eq__(self, other: object) -> bool:
        return isinstance(other, DiagramIR) and self._tup() == other._tup()

    def _tup(self) -> tuple:
        return (self.schema_version, self.object_type, self.family, self.variables, self.diagrams)

    def __repr__(self) -> str:
        return f"DiagramIR({self.family.family_id}, n={self.family.order}, {len(self.diagrams)} diagrams)"

    def to_dict(self) -> dict:
        return {
            "schema_version": self.schema_version,
            "ir_kind": "DiagramIR",
            "object_type": self.object_type,
            "family": self.family.to_dict(),
            "variables": self.variables.to_dict(),
            "diagrams": [d.to_dict() for d in self.diagrams],
        }

    @classmethod
    def from_dict(cls, d: dict) -> DiagramIR:
        import json

        from .canonical import validate_ir

        if not isinstance(d, dict) or d.get("ir_kind") != "DiagramIR":
            raise ValueError("not a DiagramIR payload")
        try:
            ir = cls(
                Family.from_dict(d["family"]),
                Variables.from_dict(d["variables"]),
                tuple(Diagram.from_dict(x) for x in d["diagrams"]),
                schema_version=d["schema_version"],
                object_type=d["object_type"],
            )
            if json.dumps(d, sort_keys=True, allow_nan=False) != json.dumps(ir.to_dict(), sort_keys=True):
                raise ValueError("unknown fields or noncanonical field types")
            validate_ir(ir)
        except (KeyError, TypeError, IndexError) as exc:
            raise ValueError("malformed DiagramIR payload") from exc
        return ir

    def to_json(self) -> str:
        import json

        return json.dumps(self.to_dict(), indent=2, sort_keys=True)

    @classmethod
    def from_json(cls, text: str) -> DiagramIR:
        import json

        return cls.from_dict(json.loads(text))
