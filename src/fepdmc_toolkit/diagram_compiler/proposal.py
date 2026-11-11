"""Minimal reversible proposal contract. Not a Monte Carlo framework."""

from __future__ import annotations

import json
import math
from dataclasses import dataclass

import numpy as np

SCHEMA_VERSION = 1
KIND_SYMMETRIC = "symmetric"
KIND_PROVIDED_LOG_RATIO = "provided_log_ratio"


class ProposalError(ValueError):
    """Illegal proposal ratio or spec."""


@dataclass(frozen=True)
class ProposalSpec:
    schema_version: int
    kind: str

    def log_q_reverse_minus_forward(self, provided: float) -> float:
        if not isinstance(provided, (int, float, np.floating)) or isinstance(provided, bool):
            raise ProposalError("log proposal ratio must be a real float")
        value = float(provided)
        if not math.isfinite(value):
            raise ProposalError("nonfinite log proposal ratio")
        if self.kind == KIND_SYMMETRIC:
            if value != 0.0:
                raise ProposalError("symmetric proposal requires log q ratio 0")
            return 0.0
        if self.kind == KIND_PROVIDED_LOG_RATIO:
            return value
        raise ProposalError(f"unknown proposal kind {self.kind}")

    def reverse_log_ratio(self, provided: float) -> float:
        return -self.log_q_reverse_minus_forward(provided)

    def to_dict(self) -> dict:
        return {"kind": self.kind, "schema_version": self.schema_version}

    def to_json(self) -> str:
        return json.dumps(self.to_dict(), indent=2, sort_keys=True) + "\n"

    @classmethod
    def from_dict(cls, payload: dict) -> ProposalSpec:
        if not isinstance(payload, dict):
            raise ProposalError("proposal payload must be an object")
        extra = set(payload) - {"kind", "schema_version"}
        if extra:
            raise ProposalError(f"unknown fields {sorted(extra)}")
        spec = cls(schema_version=int(payload["schema_version"]), kind=str(payload["kind"]))
        if spec.schema_version != SCHEMA_VERSION:
            raise ProposalError("unsupported proposal schema")
        if spec.kind not in (KIND_SYMMETRIC, KIND_PROVIDED_LOG_RATIO):
            raise ProposalError(f"unknown proposal kind {spec.kind}")
        if json.dumps(payload, sort_keys=True, allow_nan=False) != json.dumps(spec.to_dict(), sort_keys=True):
            raise ProposalError("unknown fields or noncanonical field types")
        return spec

    @classmethod
    def from_json(cls, text: str) -> ProposalSpec:
        return cls.from_dict(json.loads(text))

    @classmethod
    def symmetric(cls) -> ProposalSpec:
        return cls(SCHEMA_VERSION, KIND_SYMMETRIC)

    @classmethod
    def provided_log_ratio(cls) -> ProposalSpec:
        return cls(SCHEMA_VERSION, KIND_PROVIDED_LOG_RATIO)
