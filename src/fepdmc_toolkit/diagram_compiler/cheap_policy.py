"""Versioned cheap-approximation policy over DiagramIR categories.

Physical DiagramIR is not mutated. The policy only says which semantic
categories a cheap lowering may KEEP, DROP, APPROXIMATE, or REFUSE.
Unlisted categories are uncovered and fail closed.
"""

from __future__ import annotations

import json
from dataclasses import dataclass

POLICY_PROPAGATOR_ONLY_V1 = "propagator_only_v1"
SCHEMA_VERSION = 1

KNOWN_CATEGORIES = (
    "factor.electron_propagator",
    "factor.phonon_propagator",
    "factor.vertex",
    "factor.prefactor",
    "diagram.sign",
    "diagram.multiplicity",
    "family.global_factors",
    "diagram_sum.all_pairings",
    "momentum.q_terms",
    "electron_interface.matrix2",
    "eval.matmul",
    "eval.trace",
)

_ACTIONS = ("keep", "drop", "approximate", "refuse")


class CheapPolicyError(ValueError):
    """Unsupported, uncovered, or malformed cheap policy."""


def _tuple_of_str(values, *, field: str) -> tuple[str, ...]:
    if not isinstance(values, (list, tuple)):
        raise CheapPolicyError(f"{field} must be a list of category names")
    out = []
    for item in values:
        if not isinstance(item, str) or not item:
            raise CheapPolicyError(f"{field} entries must be nonempty strings")
        out.append(item)
    return tuple(out)


@dataclass(frozen=True)
class CheapPolicy:
    schema_version: int
    name: str
    keep: tuple[str, ...]
    drop: tuple[str, ...]
    approximate: tuple[str, ...]
    refuse: tuple[str, ...]
    notes: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        if self.schema_version != SCHEMA_VERSION:
            raise CheapPolicyError(f"unsupported policy schema {self.schema_version}")
        if not isinstance(self.name, str) or not self.name:
            raise CheapPolicyError("policy name must be a nonempty string")
        listed = self.keep + self.drop + self.approximate + self.refuse
        if len(listed) != len(set(listed)):
            raise CheapPolicyError("duplicate category in KEEP/DROP/APPROXIMATE/REFUSE")
        for category in listed:
            if category not in KNOWN_CATEGORIES:
                raise CheapPolicyError(f"unknown category {category}")

    def action(self, category: str) -> str:
        if category in self.keep:
            return "keep"
        if category in self.drop:
            return "drop"
        if category in self.approximate:
            return "approximate"
        if category in self.refuse:
            return "refuse"
        raise CheapPolicyError(f"uncovered category {category}")

    def require_allowed(self, category: str) -> str:
        act = self.action(category)
        if act == "refuse":
            raise CheapPolicyError(f"refuse category {category}")
        return act

    def present_categories(self, ir) -> tuple[str, ...]:
        from .ir import KIND_ELECTRON, MODEL_TWOBAND

        present: list[str] = [
            "diagram.sign",
            "diagram.multiplicity",
            "family.global_factors",
            "diagram_sum.all_pairings",
        ]
        seen_kind: set[str] = set()
        for diagram in ir.diagrams:
            for factor in diagram.factors:
                seen_kind.add(factor.kind)
                if factor.kind == KIND_ELECTRON:
                    present.append("momentum.q_terms")
        for kind in sorted(seen_kind):
            present.append(f"factor.{kind}")
        if ir.family.model == MODEL_TWOBAND:
            present.append("electron_interface.matrix2")
        # Cheap lowering never emits these; the policy must still declare them.
        present.extend(("eval.matmul", "eval.trace"))
        # Preserve order, drop duplicates.
        out: list[str] = []
        for category in present:
            if category not in out:
                out.append(category)
        return tuple(out)

    def check_ir(self, ir) -> None:
        for category in self.present_categories(ir):
            self.require_allowed(category)

    def to_dict(self) -> dict:
        return {
            "approximate": list(self.approximate),
            "drop": list(self.drop),
            "keep": list(self.keep),
            "name": self.name,
            "notes": list(self.notes),
            "refuse": list(self.refuse),
            "schema_version": self.schema_version,
        }

    def to_json(self) -> str:
        return json.dumps(self.to_dict(), indent=2, sort_keys=True) + "\n"

    @classmethod
    def from_dict(cls, payload: dict) -> CheapPolicy:
        if not isinstance(payload, dict):
            raise CheapPolicyError("policy payload must be an object")
        allowed = {"approximate", "drop", "keep", "name", "notes", "refuse", "schema_version"}
        extra = set(payload) - allowed
        if extra:
            raise CheapPolicyError(f"unknown fields {sorted(extra)}")
        try:
            policy = cls(
                schema_version=int(payload["schema_version"]),
                name=str(payload["name"]),
                keep=_tuple_of_str(payload.get("keep", ()), field="keep"),
                drop=_tuple_of_str(payload.get("drop", ()), field="drop"),
                approximate=_tuple_of_str(payload.get("approximate", ()), field="approximate"),
                refuse=_tuple_of_str(payload.get("refuse", ()), field="refuse"),
                notes=_tuple_of_str(payload.get("notes", ()), field="notes"),
            )
        except (KeyError, TypeError, ValueError) as exc:
            raise CheapPolicyError("malformed CheapPolicy payload") from exc
        canonical = json.dumps(payload, sort_keys=True, allow_nan=False)
        if canonical != json.dumps(policy.to_dict(), sort_keys=True):
            raise CheapPolicyError("unknown fields or noncanonical field types")
        return policy

    @classmethod
    def from_json(cls, text: str) -> CheapPolicy:
        return cls.from_dict(json.loads(text))


def propagator_only_v1() -> CheapPolicy:
    """Chosen stage-2 policy: cheap state weight from G, D, and prefactor.

    Two-band matrix electrons are replaced by scalar Gel of the same
    k_form. Vertices and matrix contraction are dropped. On the scalar
    family this policy is identical to the exact evaluator.
    """
    return CheapPolicy(
        schema_version=SCHEMA_VERSION,
        name=POLICY_PROPAGATOR_ONLY_V1,
        keep=(
            "factor.electron_propagator",
            "factor.phonon_propagator",
            "factor.prefactor",
            "diagram.sign",
            "diagram.multiplicity",
            "family.global_factors",
            "diagram_sum.all_pairings",
            "momentum.q_terms",
        ),
        drop=(
            "factor.vertex",
            "eval.matmul",
            "eval.trace",
        ),
        approximate=("electron_interface.matrix2",),
        refuse=(),
        notes=(
            "CHEAP_EVALUATOR_SEMANTICS=STATE_WEIGHT",
            "Two-band matrix G is replaced by scalar Gel of the same k_form.",
            "Empty two-band global_factors are approximated by (g/sqrt(L))**(2n).",
            "On scalar families this policy equals the exact evaluator.",
            "Historical P1 prop score is a different transition object.",
        ),
    )


def resolve_policy(policy: CheapPolicy | str | None) -> CheapPolicy:
    if policy is None or policy == POLICY_PROPAGATOR_ONLY_V1:
        return propagator_only_v1()
    if isinstance(policy, CheapPolicy):
        return policy
    raise CheapPolicyError(f"unknown policy {policy!r}")
