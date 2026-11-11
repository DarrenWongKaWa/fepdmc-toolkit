"""Independent NativeKernelIR interpreter. Does not call evaluate_transition."""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from .evaluator import Binding
from .native_ir import (
    ERROR_NAMES,
    INVALID_CHEAP_WEIGHT,
    INVALID_EXACT_TARGET,
    INVALID_IR,
    INVALID_PROPOSAL_RATIO,
    INVALID_RANDOM_UNIFORM,
    STATUS_OK,
    NativeIRError,
    NativeKernelIR,
)
from .native_validate import validate_native_ir
from .proposal import KIND_SYMMETRIC

_IMAG_ATOL = 1e-12
_RNG_FLOOR = 1e-300


@dataclass(frozen=True)
class NativeResult:
    accepted: bool
    reject_stage: int
    ell_hat: float | None
    ell_R: float | None
    exact_y_evaluated: bool
    status: int
    stage1_pass: bool
    stage2_pass: bool | None

    def to_dict(self) -> dict:
        return {
            "accepted": self.accepted,
            "ell_R": self.ell_R,
            "ell_hat": self.ell_hat,
            "exact_y_evaluated": self.exact_y_evaluated,
            "reject_stage": self.reject_stage,
            "stage1_pass": self.stage1_pass,
            "stage2_pass": self.stage2_pass,
            "status": self.status,
            "status_name": ERROR_NAMES[self.status],
        }


def _xi(k_eff: float, t: float) -> float:
    return 2.0 * t * (1.0 - math.cos(k_eff))


def _load_field(binding: Binding, field: str, index=None):
    if field == "k":
        return float(binding.k)
    if field == "t":
        return float(binding.t)
    if field == "omega":
        return float(binding.omega)
    if field == "g":
        return float(binding.g)
    if field == "L":
        return float(binding.L)
    if field == "delta":
        return float(binding.delta)
    if field == "gap":
        return float(binding.gap)
    if field == "q":
        return float(binding.q[int(index)])
    if field == "tau":
        return float(binding.tau[int(index)])
    raise NativeIRError(INVALID_IR, f"unknown field {field}")


def _as_c128(value):
    return np.complex128(value)


def _as_m2(value):
    array = np.asarray(value, dtype=np.complex128)
    if array.shape != (2, 2):
        raise NativeIRError(INVALID_IR, "expected 2x2 matrix")
    return array


def _eval_op(op, env: dict, *, binding: Binding | None, inputs: dict):
    args = [env[name] if name in env else inputs.get(name) for name in op.args]
    name = op.op
    if name == "const_f64":
        return float(op.attrs["value"])
    if name == "const_c128":
        return np.complex128(complex(op.attrs["re"], op.attrs["im"]))
    if name == "const_c128_m2":
        return np.array([[complex(re, im) for re, im in row] for row in op.attrs["value"]], dtype=np.complex128)
    if name == "load":
        return _load_field(binding, op.attrs["field"])
    if name == "load_index":
        return _load_field(binding, op.attrs["field"], op.attrs["index"])
    if name == "load_input":
        return inputs[op.attrs["name"]]
    if name == "add_f64":
        return float(args[0] + args[1])
    if name == "sub_f64":
        return float(args[0] - args[1])
    if name == "mul_f64":
        return float(args[0] * args[1])
    if name == "min_f64":
        return float(min(args[0], args[1]))
    if name == "lt_f64":
        return bool(args[0] < args[1])
    if name == "add_c128":
        return _as_c128(args[0] + args[1])
    if name == "mul_c128":
        return _as_c128(args[0] * args[1])
    if name == "add_m2":
        return _as_m2(args[0]) + _as_m2(args[1])
    if name == "scale_m2":
        return _as_m2(args[0]) * _as_c128(args[1])
    if name == "matmul":
        return _as_m2(args[0]) @ _as_m2(args[1])
    if name == "trace":
        return np.complex128(np.trace(_as_m2(args[0])))
    if name == "prim.electron_scalar":
        k_eff, t, dtau = (float(x) for x in args)
        return np.complex128(math.exp(-_xi(k_eff, t) * dtau))
    if name == "prim.electron_twoband":
        k_eff, t, dtau, delta, gap = (float(x) for x in args)
        xi = _xi(k_eff, t)
        h = np.array([[xi, delta], [delta, xi + gap]], dtype=np.complex128)
        w, v = np.linalg.eigh(h)
        return (v * np.exp(-w * dtau)) @ v.conj().T
    if name == "prim.phonon":
        omega, dtau = float(args[0]), float(args[1])
        return np.complex128(math.exp(-omega * dtau))
    if name == "prim.prefactor":
        g, L = float(args[0]), float(args[1])
        n = int(op.attrs["n"])
        return np.complex128((g / math.sqrt(L)) ** (2 * n))
    if name == "prim.vertex_sigmaz":
        g, L = float(args[0]), float(args[1])
        scale = g / math.sqrt(L)
        return scale * np.array([[1.0, 0.0], [0.0, -1.0]], dtype=np.complex128)
    if name == "log_positive_real":
        number = _as_c128(args[0])
        code = INVALID_CHEAP_WEIGHT if op.attrs.get("on_fail") == "INVALID_CHEAP_WEIGHT" else INVALID_EXACT_TARGET
        imag_atol = float(op.attrs.get("imag_atol", _IMAG_ATOL))
        if not np.isfinite(number.real) or not np.isfinite(number.imag):
            raise NativeIRError(code, "nonfinite weight")
        if abs(number.imag) > imag_atol or number.real <= 0.0:
            raise NativeIRError(code, "non-positive or complex weight")
        return float(math.log(number.real))
    if name == "log_uniform":
        which = op.attrs["which"]
        u = float(inputs[which])
        floor = float(op.attrs.get("floor", _RNG_FLOOR))
        if isinstance(inputs[which], bool) or not math.isfinite(u) or u < 0.0 or u > 1.0:
            raise NativeIRError(INVALID_RANDOM_UNIFORM, "uniform not in [0,1]")
        return float(math.log(max(u, floor)))
    if name == "check_log_q":
        value = float(inputs["log_q"])
        if not math.isfinite(value):
            raise NativeIRError(INVALID_PROPOSAL_RATIO, "nonfinite log q")
        if op.attrs.get("kind") == KIND_SYMMETRIC and value != 0.0:
            raise NativeIRError(INVALID_PROPOSAL_RATIO, "symmetric requires 0")
        return value
    raise NativeIRError(INVALID_IR, f"unsupported op {name}")


def _run_graph(graph, binding: Binding):
    env = {}
    for op in graph.ops:
        local = dict(op.attrs)
        if local.get("binding") == "binding":
            # graph parameter; actual Binding already selected
            pass
        env[op.id] = _eval_op(op, env, binding=binding, inputs={})
    return env[graph.output]


def interpret_native_kernel(
    nir: NativeKernelIR,
    *,
    x: Binding | None = None,
    y: Binding | None = None,
    log_q: float = 0.0,
    u1: float = 0.5,
    u2: float = 0.5,
    exact_x_valid: bool = False,
    exact_log_weight_x: float = 0.0,
    log_wx=None,
    log_wy=None,
    log_px=None,
    log_py=None,
) -> NativeResult:
    validate_native_ir(nir)
    try:
        return _interpret_body(
            nir,
            x=x,
            y=y,
            log_q=log_q,
            u1=u1,
            u2=u2,
            exact_x_valid=exact_x_valid,
            exact_log_weight_x=exact_log_weight_x,
            log_wx=log_wx,
            log_wy=log_wy,
            log_px=log_px,
            log_py=log_py,
        )
    except NativeIRError as exc:
        if exc.code == INVALID_IR:
            raise
        return NativeResult(
            accepted=False,
            reject_stage=0,
            ell_hat=None,
            ell_R=None,
            exact_y_evaluated=False,
            status=exc.code,
            stage1_pass=False,
            stage2_pass=None,
        )


def _interpret_body(
    nir: NativeKernelIR,
    *,
    x,
    y,
    log_q,
    u1,
    u2,
    exact_x_valid,
    exact_log_weight_x,
    log_wx,
    log_wy,
    log_px,
    log_py,
) -> NativeResult:
    inputs = {
        "in.exact_x_valid": bool(exact_x_valid),
        "in.exact_log_weight_x": float(exact_log_weight_x),
        "in.log_q": float(log_q),
        "u1": float(u1),
        "u2": float(u2),
        "log_q": float(log_q),
        "exact_log_weight_x": float(exact_log_weight_x),
        "log_wx": log_wx,
        "log_wy": log_wy,
        "log_px": log_px,
        "log_py": log_py,
    }
    env: dict = {"in.exact_x_valid": bool(exact_x_valid)}
    block_name = nir.entry
    pred = None
    ell_hat = None
    ell_R = None
    while True:
        block = nir.block(block_name)
        for op in block.ops:
            if op.op == "call_graph":
                graph = nir.graph(op.attrs["graph"])
                side = op.attrs["binding"]
                binding = x if side == "x" else y
                env[op.id] = _run_graph(graph, binding)
            elif op.op == "phi":
                chosen = None
                for src_block, src_id in op.attrs["preds"]:
                    if src_block == pred:
                        chosen = env[src_id]
                if chosen is None:
                    raise NativeIRError(INVALID_IR, "phi has no predecessor")
                env[op.id] = chosen
            else:
                env[op.id] = _eval_op(op, env, binding=None, inputs=inputs)
            if op.id == "ell_hat":
                ell_hat = float(env[op.id])
            if op.id == "ell_R":
                ell_R = float(env[op.id])
        term = block.term
        if term.op == "br":
            pred, block_name = block_name, str(term.attrs["target"])
            continue
        if term.op == "br_cond":
            cond_name = term.args[0] if term.args else None
            if cond_name == "in.exact_x_valid":
                flag = bool(exact_x_valid)
            else:
                flag = bool(env[cond_name])
            pred, block_name = block_name, str(term.attrs["true" if flag else "false"])
            continue
        if term.op == "return":
            accepted = bool(term.attrs["accepted"])
            reject_stage = int(term.attrs["reject_stage"])
            exact_y = bool(term.attrs["exact_y_evaluated"])
            ell_R_valid = bool(term.attrs["ell_R_valid"])
            return NativeResult(
                accepted=accepted,
                reject_stage=reject_stage,
                ell_hat=ell_hat,
                ell_R=ell_R if ell_R_valid else None,
                exact_y_evaluated=exact_y,
                status=STATUS_OK,
                stage1_pass=reject_stage != 1,
                stage2_pass=True if accepted else (False if reject_stage == 2 else None),
            )
        raise NativeIRError(INVALID_IR, f"unknown terminator {term.op}")


def native_acceptance_probability(log_wx, log_wy, log_px, log_py, log_q) -> float:
    """Execute Design B score ops (same formulas as NativeKernelIR CFG)."""
    ell_hat = float(log_wy - log_wx)
    ell_R = float((log_py - log_px) + log_q)
    a1 = 1.0 if ell_hat >= 0.0 else math.exp(ell_hat)
    delta = ell_R - ell_hat
    a2 = 1.0 if delta >= 0.0 else math.exp(delta)
    return float(a1 * a2)
