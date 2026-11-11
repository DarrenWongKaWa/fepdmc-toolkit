"""Versioned exact-target policy for delayed acceptance.

pi(x) is the stage-1 exact F when F is finite, real, and strictly positive.
|F|, |F|_C, and F^2 are refused. This is a restricted positive-real subset,
not a new Monte Carlo measure.
"""

from __future__ import annotations

import json
import math
from dataclasses import dataclass

import numpy as np

SCHEMA_VERSION = 1
POLICY_POSITIVE_REAL_F_V1 = "positive_real_F_v1"
_IMAG_ATOL = 1e-12
_REFUSE_TRANSFORMS = ("abs_F", "F_squared", "complex_modulus")


class TargetPolicyError(ValueError):
    """Exact F cannot be used as a positive target weight."""


@dataclass(frozen=True)
class TargetPolicy:
    schema_version: int
    name: str
    target: str
    refuse_transforms: tuple[str, ...]
    imag_atol: float
    notes: tuple[str, ...] = ()

    def log_weight(self, value) -> float:
        if self.target != "exact_F":
            raise TargetPolicyError(f"unsupported target {self.target}")
        number = complex(value)
        if not np.isfinite(number.real) or not np.isfinite(number.imag):
            raise TargetPolicyError("nonfinite exact F")
        if abs(number.imag) > self.imag_atol:
            raise TargetPolicyError("complex exact F")
        if number.real <= 0.0:
            raise TargetPolicyError("non-positive exact F")
        return float(math.log(number.real))

    def to_dict(self) -> dict:
        return {
            "imag_atol": self.imag_atol,
            "name": self.name,
            "notes": list(self.notes),
            "refuse_transforms": list(self.refuse_transforms),
            "schema_version": self.schema_version,
            "target": self.target,
        }

    def to_json(self) -> str:
        return json.dumps(self.to_dict(), indent=2, sort_keys=True) + "\n"

    @classmethod
    def from_dict(cls, payload: dict) -> TargetPolicy:
        if not isinstance(payload, dict):
            raise TargetPolicyError("target policy payload must be an object")
        allowed = {"imag_atol", "name", "notes", "refuse_transforms", "schema_version", "target"}
        extra = set(payload) - allowed
        if extra:
            raise TargetPolicyError(f"unknown fields {sorted(extra)}")
        policy = cls(
            schema_version=int(payload["schema_version"]),
            name=str(payload["name"]),
            target=str(payload["target"]),
            refuse_transforms=tuple(payload.get("refuse_transforms", ())),
            imag_atol=float(payload["imag_atol"]),
            notes=tuple(payload.get("notes", ())),
        )
        if json.dumps(payload, sort_keys=True, allow_nan=False) != json.dumps(policy.to_dict(), sort_keys=True):
            raise TargetPolicyError("unknown fields or noncanonical field types")
        if policy.schema_version != SCHEMA_VERSION:
            raise TargetPolicyError("unsupported target-policy schema")
        return policy

    @classmethod
    def from_json(cls, text: str) -> TargetPolicy:
        return cls.from_dict(json.loads(text))


def positive_real_F_v1() -> TargetPolicy:
    return TargetPolicy(
        schema_version=SCHEMA_VERSION,
        name=POLICY_POSITIVE_REAL_F_V1,
        target="exact_F",
        refuse_transforms=_REFUSE_TRANSFORMS,
        imag_atol=_IMAG_ATOL,
        notes=(
            "pi(x)=F(x) on the restricted set F real, finite, and strictly positive.",
            "Scalar Holstein g!=0 is in this set. Two-band n=1 is positive for g!=0; n>=2 is admitted only when the sample F>0.",
            "Do not sample |F| or F^2.",
        ),
    )
