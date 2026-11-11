"""Design B delayed-acceptance kernel over stage-1 exact and stage-2 cheap evaluators.

ell_hat = log W_hat(y) - log W_hat(x)          # stage-2 STATE_WEIGHT; no q
ell_R   = log pi(y) - log pi(x) + log q(y,x) - log q(x,y)

Stage 1 uses ell_hat only. Stage 2 repairs both target approximation and
Hastings. Not a general Monte Carlo framework. No Fortran, no P1.
"""

from __future__ import annotations

import hashlib
import json
import math
from dataclasses import dataclass

import numpy as np

from .canonical import validate_ir
from .cheap_evaluator import CheapEvaluator, compile_cheap, compile_exact
from .evaluator import Binding, validate_binding
from .ir import DiagramIR
from .proposal import ProposalError, ProposalSpec
from .target_policy import TargetPolicy, positive_real_F_v1

SCHEMA_VERSION = 1
DESIGN_B = "B"
RNG_LOG_FLOOR = 1e-300
ATOL = 1e-12
RTOL = 1e-10


class DAError(ValueError):
    """Illegal delayed-acceptance query or specification."""


@dataclass(frozen=True)
class NumericalPolicy:
    rng_log_floor: float = RNG_LOG_FLOOR
    clip: str = "none"
    atol: float = ATOL
    rtol: float = RTOL

    def to_dict(self) -> dict:
        return {
            "atol": self.atol,
            "clip": self.clip,
            "rng_log_floor": self.rng_log_floor,
            "rtol": self.rtol,
        }


@dataclass(frozen=True)
class DelayedAcceptanceKernelSpec:
    schema_version: int
    spec_kind: str
    design: str
    diagram_ir_schema: str
    diagram_ir_digest: str
    exact_lowering: str
    cheap_policy_name: str
    cheap_policy_schema: int
    target_policy_name: str
    proposal_kind: str
    numerical: dict
    stage1: str
    stage2: str

    def to_dict(self) -> dict:
        return {
            "cheap_policy_name": self.cheap_policy_name,
            "cheap_policy_schema": self.cheap_policy_schema,
            "design": self.design,
            "diagram_ir_digest": self.diagram_ir_digest,
            "diagram_ir_schema": self.diagram_ir_schema,
            "exact_lowering": self.exact_lowering,
            "numerical": dict(self.numerical),
            "proposal_kind": self.proposal_kind,
            "schema_version": self.schema_version,
            "spec_kind": self.spec_kind,
            "stage1": self.stage1,
            "stage2": self.stage2,
            "target_policy_name": self.target_policy_name,
        }

    def to_json(self) -> str:
        return json.dumps(self.to_dict(), indent=2, sort_keys=True) + "\n"

    @classmethod
    def from_dict(cls, payload: dict) -> DelayedAcceptanceKernelSpec:
        if not isinstance(payload, dict) or payload.get("spec_kind") != "DelayedAcceptanceKernelSpec":
            raise DAError("not a DelayedAcceptanceKernelSpec payload")
        allowed = set(
            cls(
                schema_version=1,
                spec_kind="DelayedAcceptanceKernelSpec",
                design="B",
                diagram_ir_schema="",
                diagram_ir_digest="",
                exact_lowering="",
                cheap_policy_name="",
                cheap_policy_schema=1,
                target_policy_name="",
                proposal_kind="",
                numerical={},
                stage1="",
                stage2="",
            ).to_dict()
        )
        extra = set(payload) - allowed
        if extra:
            raise DAError(f"unknown fields {sorted(extra)}")
        spec = cls(
            schema_version=int(payload["schema_version"]),
            spec_kind=str(payload["spec_kind"]),
            design=str(payload["design"]),
            diagram_ir_schema=str(payload["diagram_ir_schema"]),
            diagram_ir_digest=str(payload["diagram_ir_digest"]),
            exact_lowering=str(payload["exact_lowering"]),
            cheap_policy_name=str(payload["cheap_policy_name"]),
            cheap_policy_schema=int(payload["cheap_policy_schema"]),
            target_policy_name=str(payload["target_policy_name"]),
            proposal_kind=str(payload["proposal_kind"]),
            numerical=dict(payload["numerical"]),
            stage1=str(payload["stage1"]),
            stage2=str(payload["stage2"]),
        )
        if spec.schema_version != SCHEMA_VERSION or spec.design != DESIGN_B:
            raise DAError("unsupported DA kernel spec")
        if json.dumps(payload, sort_keys=True, allow_nan=False) != json.dumps(spec.to_dict(), sort_keys=True):
            raise DAError("unknown fields or noncanonical field types")
        return spec

    @classmethod
    def from_json(cls, text: str) -> DelayedAcceptanceKernelSpec:
        return cls.from_dict(json.loads(text))


@dataclass(frozen=True)
class DATransitionResult:
    cheap_log_weight_x: float
    cheap_log_weight_y: float
    ell_hat: float
    exact_log_weight_x: float | None
    exact_log_weight_y: float | None
    ell_R: float | None
    log_q_reverse_minus_forward: float
    stage1_pass: bool
    exact_y_evaluated: bool
    exact_x_evaluated: bool
    stage2_pass: bool | None
    accepted: bool


def ell_hat_state_weight(log_w_x: float, log_w_y: float) -> float:
    return float(log_w_y - log_w_x)


def ell_R_design_b(log_pi_x: float, log_pi_y: float, log_q_reverse_minus_forward: float) -> float:
    return float((log_pi_y - log_pi_x) + log_q_reverse_minus_forward)


def acceptance_probability(ell_hat: float, ell_R: float) -> float:
    alpha1 = 1.0 if ell_hat >= 0.0 else math.exp(ell_hat)
    delta = ell_R - ell_hat
    alpha2 = 1.0 if delta >= 0.0 else math.exp(delta)
    return float(alpha1 * alpha2)


def delayed_acceptance_identity(
    a_xy: float, a_yx: float, r_exact: float, *, atol: float = ATOL, rtol: float = RTOL
) -> bool:
    return bool(abs(a_xy - r_exact * a_yx) <= atol + rtol * abs(r_exact * a_yx))


def log_uniform(u: float, floor: float = RNG_LOG_FLOOR) -> float:
    if isinstance(u, bool) or not isinstance(u, (int, float, np.floating)):
        raise DAError("uniform must be a real float")
    value = float(u)
    if not math.isfinite(value) or value < 0.0 or value > 1.0:
        raise DAError("uniform must lie in [0, 1]")
    return float(math.log(max(value, floor)))


def stage1_pass(ell_hat: float, u1: float, floor: float = RNG_LOG_FLOOR) -> bool:
    return log_uniform(u1, floor) < min(0.0, ell_hat)


def stage2_pass(ell_R: float, ell_hat: float, u2: float, floor: float = RNG_LOG_FLOOR) -> bool:
    return log_uniform(u2, floor) < min(0.0, ell_R - ell_hat)


def transition_matrix(
    pi,
    cheap,
    q,
    *,
    stage2: bool = True,
    drop_ell_hat_in_stage2: bool = False,
    invert_hastings: bool = False,
):
    pi = np.asarray(pi, dtype=np.float64)
    cheap = np.asarray(cheap, dtype=np.float64)
    q = np.asarray(q, dtype=np.float64)
    n = pi.shape[0]
    accept = np.zeros((n, n), dtype=np.float64)
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            ell_hat = ell_hat_state_weight(math.log(cheap[i]), math.log(cheap[j]))
            log_q = math.log(q[j, i]) - math.log(q[i, j])
            if invert_hastings:
                log_q = -log_q
            ell_R = ell_R_design_b(math.log(pi[i]), math.log(pi[j]), log_q)
            if not stage2:
                accept[i, j] = 1.0 if ell_hat >= 0.0 else math.exp(ell_hat)
            elif drop_ell_hat_in_stage2:
                a1 = 1.0 if ell_hat >= 0.0 else math.exp(ell_hat)
                a2 = 1.0 if ell_R >= 0.0 else math.exp(ell_R)
                accept[i, j] = a1 * a2
            else:
                accept[i, j] = acceptance_probability(ell_hat, ell_R)
    p = q * accept
    for i in range(n):
        p[i, i] = 1.0 - float(np.sum(p[i]))
    return p


def _ir_digest(ir: DiagramIR) -> str:
    payload = json.dumps(ir.to_dict(), sort_keys=True, separators=(",", ":"), allow_nan=False)
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


class DelayedAcceptanceKernel:
    def __init__(
        self,
        ir: DiagramIR,
        exact,
        cheap: CheapEvaluator,
        target_policy: TargetPolicy,
        proposal_spec: ProposalSpec,
        numerical: NumericalPolicy,
        spec: DelayedAcceptanceKernelSpec,
    ) -> None:
        self.ir = ir
        self._exact = exact
        self._cheap = cheap
        self._target = target_policy
        self._proposal = proposal_spec
        self._numerical = numerical
        self.spec = spec

    def cheap_ell_hat(self, current: Binding, proposed: Binding) -> float:
        return float(self._cheap.transition_score(current, proposed))

    def exact_log_weight(self, binding: Binding) -> float:
        validate_binding(self.ir, binding)
        value = self._exact.evaluate(binding)["F"]
        return self._target.log_weight(value)

    def evaluate_transition(
        self,
        current: Binding,
        proposed: Binding,
        log_q_reverse_minus_forward: float,
        u1: float,
        u2: float,
        exact_log_weight_x: float | None = None,
        exact_x: Binding | None = None,
    ) -> DATransitionResult:
        validate_binding(self.ir, current)
        validate_binding(self.ir, proposed)
        try:
            log_q = self._proposal.log_q_reverse_minus_forward(log_q_reverse_minus_forward)
        except ProposalError as exc:
            raise DAError(str(exc)) from exc
        log_wx = self._cheap.log_weight(current)
        log_wy = self._cheap.log_weight(proposed)
        ell_hat = ell_hat_state_weight(log_wx, log_wy)
        if not math.isfinite(ell_hat):
            raise DAError("nonfinite ell_hat")
        s1 = stage1_pass(ell_hat, u1, self._numerical.rng_log_floor)
        if not s1:
            return DATransitionResult(
                cheap_log_weight_x=log_wx,
                cheap_log_weight_y=log_wy,
                ell_hat=ell_hat,
                exact_log_weight_x=None,
                exact_log_weight_y=None,
                ell_R=None,
                log_q_reverse_minus_forward=log_q,
                stage1_pass=False,
                exact_y_evaluated=False,
                exact_x_evaluated=False,
                stage2_pass=None,
                accepted=False,
            )
        exact_x_evaluated = False
        if exact_log_weight_x is not None:
            if exact_x is None or exact_x != current:
                raise DAError("stale or unkeyed exact(x) cache")
            if not math.isfinite(float(exact_log_weight_x)):
                raise DAError("nonfinite cached exact log-weight")
            log_px = float(exact_log_weight_x)
        else:
            log_px = self.exact_log_weight(current)
            exact_x_evaluated = True
        log_py = self.exact_log_weight(proposed)
        ell_R = ell_R_design_b(log_px, log_py, log_q)
        if not math.isfinite(ell_R):
            raise DAError("nonfinite ell_R")
        s2 = stage2_pass(ell_R, ell_hat, u2, self._numerical.rng_log_floor)
        return DATransitionResult(
            cheap_log_weight_x=log_wx,
            cheap_log_weight_y=log_wy,
            ell_hat=ell_hat,
            exact_log_weight_x=log_px,
            exact_log_weight_y=log_py,
            ell_R=ell_R,
            log_q_reverse_minus_forward=log_q,
            stage1_pass=True,
            exact_y_evaluated=True,
            exact_x_evaluated=exact_x_evaluated,
            stage2_pass=s2,
            accepted=bool(s2),
        )


def compile_delayed_acceptance(
    ir: DiagramIR,
    exact=None,
    cheap: CheapEvaluator | None = None,
    target_policy: TargetPolicy | None = None,
    proposal_spec: ProposalSpec | None = None,
    numerical_policy: NumericalPolicy | None = None,
) -> DelayedAcceptanceKernel:
    validate_ir(ir)
    snapshot = DiagramIR.from_json(ir.to_json())
    resolved_exact = exact if exact is not None else compile_exact(snapshot)
    resolved_cheap = cheap if cheap is not None else compile_cheap(snapshot)
    resolved_target = target_policy if target_policy is not None else positive_real_F_v1()
    resolved_proposal = proposal_spec if proposal_spec is not None else ProposalSpec.provided_log_ratio()
    resolved_numerical = numerical_policy if numerical_policy is not None else NumericalPolicy()
    if resolved_numerical.clip != "none":
        raise DAError("clipping is not part of the stage-3 default kernel")
    if resolved_cheap.policy.name != "propagator_only_v1":
        raise DAError("stage-3 kernel requires stage-2 propagator_only_v1")
    if resolved_target.name != "positive_real_F_v1" or resolved_target.target != "exact_F":
        raise DAError("stage-3 kernel requires positive_real_F_v1; |F| is refused")
    if "abs_F" not in resolved_target.refuse_transforms:
        raise DAError("stage-3 kernel refuses abs_F targets")
    spec = DelayedAcceptanceKernelSpec(
        schema_version=SCHEMA_VERSION,
        spec_kind="DelayedAcceptanceKernelSpec",
        design=DESIGN_B,
        diagram_ir_schema=snapshot.schema_version,
        diagram_ir_digest=_ir_digest(snapshot),
        exact_lowering="compile_exact/lower_diagramwise_share",
        cheap_policy_name=resolved_cheap.policy.name,
        cheap_policy_schema=resolved_cheap.policy.schema_version,
        target_policy_name=resolved_target.name,
        proposal_kind=resolved_proposal.kind,
        numerical=resolved_numerical.to_dict(),
        stage1="min(1, exp(ell_hat)); ell_hat = log W_hat(y) - log W_hat(x)",
        stage2="min(1, exp(ell_R - ell_hat)); ell_R = log pi(y)-log pi(x) + log q(y,x)-log q(x,y)",
    )
    return DelayedAcceptanceKernel(
        snapshot,
        resolved_exact,
        resolved_cheap,
        resolved_target,
        resolved_proposal,
        resolved_numerical,
        spec,
    )
