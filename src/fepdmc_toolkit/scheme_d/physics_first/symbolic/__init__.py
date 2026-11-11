"""Symbolic R1: compile-once, bind-many. Not a native D(C)."""

from .cost import NOT_COMPUTED, BuildReport, CallLedger
from .plan_spec import PREREGISTERED_SCALAR, PlanSpec, scalar_spec, twoband_spec
from .r1 import BuiltPlan, PlanCache, r1_symbolic_build
from .session import PhysicsSession, create_session

__all__ = [
    "NOT_COMPUTED",
    "BuildReport",
    "CallLedger",
    "PlanSpec",
    "PREREGISTERED_SCALAR",
    "scalar_spec",
    "twoband_spec",
    "BuiltPlan",
    "PlanCache",
    "r1_symbolic_build",
    "PhysicsSession",
    "create_session",
]
