"""Lower stage-1/2 DAGs and stage-3 Design B control flow into NativeKernelIR."""

from __future__ import annotations

import numpy as np

from .graph import lower_diagramwise
from .ir import KIND_ELECTRON, KIND_PHONON, KIND_PREFACTOR, KIND_VERTEX, MODEL_TWOBAND, MomentumForm
from .native_ir import (
    IR_KIND,
    REGION_CHEAP,
    REGION_CONTROL,
    REGION_EXACT,
    SCHEMA_VERSION,
    TYPE_BOOL,
    TYPE_C128,
    TYPE_C128_M2,
    TYPE_F64,
    BindingLayout,
    NativeBlock,
    NativeGraph,
    NativeKernelIR,
    NativeOp,
    NativeTerminator,
    binding_layout_for,
)


class _Seq:
    def __init__(self, prefix: str) -> None:
        self.prefix = prefix
        self.n = 0
        self.ops: list[NativeOp] = []

    def emit(self, op: str, type_: str, *, args=(), attrs=None, region: str) -> str:
        self.n += 1
        vid = f"{self.prefix}{self.n}"
        self.ops.append(NativeOp(vid, op, type_, tuple(args), dict(attrs or {}), region))
        return vid


def _c128_const(seq: _Seq, value, region: str) -> str:
    number = complex(value)
    return seq.emit("const_c128", TYPE_C128, attrs={"re": float(number.real), "im": float(number.imag)}, region=region)


def _m2_const(seq: _Seq, matrix, region: str) -> str:
    array = np.asarray(matrix, dtype=np.complex128)
    payload = [[[float(z.real), float(z.imag)] for z in row] for row in array]
    return seq.emit("const_c128_m2", TYPE_C128_M2, attrs={"value": payload}, region=region)


def _k_eff(seq: _Seq, binding: str, k_form: MomentumForm, region: str) -> str:
    k = seq.emit("load", TYPE_F64, attrs={"binding": binding, "field": "k"}, region=region)
    acc = seq.emit(
        "mul_f64",
        TYPE_F64,
        args=(seq.emit("const_f64", TYPE_F64, attrs={"value": float(k_form.k_coeff)}, region=region), k),
        region=region,
    )
    for slot, coeff in k_form.q_terms:
        q = seq.emit(
            "load_index", TYPE_F64, attrs={"binding": binding, "field": "q", "index": int(slot)}, region=region
        )
        term = seq.emit(
            "mul_f64",
            TYPE_F64,
            args=(seq.emit("const_f64", TYPE_F64, attrs={"value": float(coeff)}, region=region), q),
            region=region,
        )
        acc = seq.emit("add_f64", TYPE_F64, args=(acc, term), region=region)
    return acc


def _dtau(seq: _Seq, binding: str, i0: int, i1: int, region: str) -> str:
    t1 = seq.emit("load_index", TYPE_F64, attrs={"binding": binding, "field": "tau", "index": int(i1)}, region=region)
    t0 = seq.emit("load_index", TYPE_F64, attrs={"binding": binding, "field": "tau", "index": int(i0)}, region=region)
    return seq.emit("sub_f64", TYPE_F64, args=(t1, t0), region=region)


def _lower_electron(seq: _Seq, binding: str, leaf: tuple, region: str, twoband: bool) -> str:
    _, interval, _slots, k_form = leaf
    k_eff = _k_eff(seq, binding, k_form, region)
    t = seq.emit("load", TYPE_F64, attrs={"binding": binding, "field": "t"}, region=region)
    dtau = _dtau(seq, binding, interval[0], interval[1], region)
    if twoband and region == REGION_EXACT:
        delta = seq.emit("load", TYPE_F64, attrs={"binding": binding, "field": "delta"}, region=region)
        gap = seq.emit("load", TYPE_F64, attrs={"binding": binding, "field": "gap"}, region=region)
        return seq.emit("prim.electron_twoband", TYPE_C128_M2, args=(k_eff, t, dtau, delta, gap), region=region)
    return seq.emit("prim.electron_scalar", TYPE_C128, args=(k_eff, t, dtau), region=region)


def _lower_phonon(seq: _Seq, binding: str, leaf: tuple, region: str) -> str:
    _, _slot, tau_open, tau_close = leaf
    omega = seq.emit("load", TYPE_F64, attrs={"binding": binding, "field": "omega"}, region=region)
    dtau = _dtau(seq, binding, tau_open, tau_close, region)
    return seq.emit("prim.phonon", TYPE_C128, args=(omega, dtau), region=region)


def _lower_prefactor(seq: _Seq, binding: str, order: int, region: str) -> str:
    g = seq.emit("load", TYPE_F64, attrs={"binding": binding, "field": "g"}, region=region)
    L = seq.emit("load", TYPE_F64, attrs={"binding": binding, "field": "L"}, region=region)
    return seq.emit("prim.prefactor", TYPE_C128, args=(g, L), attrs={"n": int(order)}, region=region)


def _lower_vertex(seq: _Seq, binding: str, region: str) -> str:
    g = seq.emit("load", TYPE_F64, attrs={"binding": binding, "field": "g"}, region=region)
    L = seq.emit("load", TYPE_F64, attrs={"binding": binding, "field": "L"}, region=region)
    return seq.emit("prim.vertex_sigmaz", TYPE_C128_M2, args=(g, L), region=region)


def lower_cheap_graph(dag, *, prefix: str = "c") -> NativeGraph:
    seq = _Seq(prefix)
    ids: dict[int, str] = {}
    types: dict[str, str] = {}
    order = dag.order
    for index, node in enumerate(dag.nodes):
        if node.op == "const":
            vid = _c128_const(seq, node.const, REGION_CHEAP)
        elif node.op == "input":
            kind = node.leaf[0]
            if kind == KIND_ELECTRON:
                vid = _lower_electron(seq, "binding", node.leaf, REGION_CHEAP, twoband=False)
            elif kind == KIND_PHONON:
                vid = _lower_phonon(seq, "binding", node.leaf, REGION_CHEAP)
            elif kind == KIND_PREFACTOR:
                vid = _lower_prefactor(seq, "binding", order, REGION_CHEAP)
            else:
                raise ValueError(kind)
        elif node.op == "add":
            vid = seq.emit("add_c128", TYPE_C128, args=(ids[node.inputs[0]], ids[node.inputs[1]]), region=REGION_CHEAP)
        elif node.op == "mul":
            vid = seq.emit("mul_c128", TYPE_C128, args=(ids[node.inputs[0]], ids[node.inputs[1]]), region=REGION_CHEAP)
        else:
            raise ValueError(node.op)
        ids[index] = vid
        types[vid] = TYPE_C128
    return NativeGraph("cheap", REGION_CHEAP, tuple(seq.ops), ids[dag.outputs["F_hat"]], TYPE_C128)


def lower_exact_graph(dag, *, prefix: str = "e") -> NativeGraph:
    seq = _Seq(prefix)
    ids: dict[int, str] = {}
    types: dict[str, str] = {}
    twoband = dag.family_model == MODEL_TWOBAND
    order = dag.order
    for index, node in enumerate(dag.nodes):
        if node.op == "const":
            value = node.const
            if np.ndim(value) == 2:
                vid = _m2_const(seq, value, REGION_EXACT)
                types[vid] = TYPE_C128_M2
            else:
                vid = _c128_const(seq, value, REGION_EXACT)
                types[vid] = TYPE_C128
        elif node.op == "input":
            kind = node.leaf[0]
            if kind == KIND_ELECTRON:
                vid = _lower_electron(seq, "binding", node.leaf, REGION_EXACT, twoband=twoband)
            elif kind == KIND_PHONON:
                vid = _lower_phonon(seq, "binding", node.leaf, REGION_EXACT)
            elif kind == KIND_PREFACTOR:
                vid = _lower_prefactor(seq, "binding", order, REGION_EXACT)
            elif kind == KIND_VERTEX:
                vid = _lower_vertex(seq, "binding", REGION_EXACT)
            else:
                raise ValueError(kind)
            types[vid] = TYPE_C128_M2 if vid in {op.id for op in seq.ops if op.type == TYPE_C128_M2} else TYPE_C128
            types[vid] = seq.ops[-1].type
        elif node.op in ("add", "mul", "matmul"):
            a, b = ids[node.inputs[0]], ids[node.inputs[1]]
            ta, tb = types[a], types[b]
            if node.op == "matmul" or (ta == TYPE_C128_M2 and tb == TYPE_C128_M2 and node.op == "mul"):
                if node.op == "matmul" and ta == TYPE_C128_M2 and tb == TYPE_C128_M2:
                    vid = seq.emit("matmul", TYPE_C128_M2, args=(a, b), region=REGION_EXACT)
                    types[vid] = TYPE_C128_M2
                elif ta == TYPE_C128_M2 and tb == TYPE_C128:
                    vid = seq.emit("scale_m2", TYPE_C128_M2, args=(a, b), region=REGION_EXACT)
                    types[vid] = TYPE_C128_M2
                elif ta == TYPE_C128 and tb == TYPE_C128_M2:
                    vid = seq.emit("scale_m2", TYPE_C128_M2, args=(b, a), region=REGION_EXACT)
                    types[vid] = TYPE_C128_M2
                elif node.op == "add" and ta == TYPE_C128_M2:
                    vid = seq.emit("add_m2", TYPE_C128_M2, args=(a, b), region=REGION_EXACT)
                    types[vid] = TYPE_C128_M2
                else:
                    vid = seq.emit(
                        "mul_c128" if node.op != "add" else "add_c128", TYPE_C128, args=(a, b), region=REGION_EXACT
                    )
                    types[vid] = TYPE_C128
            elif node.op == "add":
                if ta == TYPE_C128_M2:
                    vid = seq.emit("add_m2", TYPE_C128_M2, args=(a, b), region=REGION_EXACT)
                    types[vid] = TYPE_C128_M2
                else:
                    vid = seq.emit("add_c128", TYPE_C128, args=(a, b), region=REGION_EXACT)
                    types[vid] = TYPE_C128
            else:
                if ta == TYPE_C128_M2 or tb == TYPE_C128_M2:
                    matrix, scalar = (a, b) if ta == TYPE_C128_M2 else (b, a)
                    vid = seq.emit("scale_m2", TYPE_C128_M2, args=(matrix, scalar), region=REGION_EXACT)
                    types[vid] = TYPE_C128_M2
                else:
                    vid = seq.emit("mul_c128", TYPE_C128, args=(a, b), region=REGION_EXACT)
                    types[vid] = TYPE_C128
        elif node.op == "trace":
            vid = seq.emit("trace", TYPE_C128, args=(ids[node.inputs[0]],), region=REGION_EXACT)
            types[vid] = TYPE_C128
        else:
            raise ValueError(node.op)
        ids[index] = vid
    out_name = "F" if "F" in dag.outputs else next(iter(dag.outputs))
    return NativeGraph("exact", REGION_EXACT, tuple(seq.ops), ids[dag.outputs[out_name]], TYPE_C128)


def _control_blocks(source: dict) -> tuple[NativeBlock, ...]:
    def op(oid, name, type_, args=(), attrs=None, region=REGION_CONTROL):
        return NativeOp(oid, name, type_, tuple(args), dict(attrs or {}), region)

    entry = NativeBlock(
        "entry",
        (
            op("cheap_x", "call_graph", TYPE_C128, attrs={"graph": "cheap", "binding": "x"}, region=REGION_CHEAP),
            op("cheap_y", "call_graph", TYPE_C128, attrs={"graph": "cheap", "binding": "y"}, region=REGION_CHEAP),
            op(
                "log_wx",
                "log_positive_real",
                TYPE_F64,
                ("cheap_x",),
                {"on_fail": "INVALID_CHEAP_WEIGHT", "imag_atol": 1e-12},
            ),
            op(
                "log_wy",
                "log_positive_real",
                TYPE_F64,
                ("cheap_y",),
                {"on_fail": "INVALID_CHEAP_WEIGHT", "imag_atol": 1e-12},
            ),
            op("ell_hat", "sub_f64", TYPE_F64, ("log_wy", "log_wx")),
            op(
                "log_q",
                "check_log_q",
                TYPE_F64,
                (),
                {"on_fail": "INVALID_PROPOSAL_RATIO", "kind": source["proposal_kind"]},
            ),
            op(
                "log_u1",
                "log_uniform",
                TYPE_F64,
                (),
                {"which": "u1", "on_fail": "INVALID_RANDOM_UNIFORM", "floor": 1e-300},
            ),
            op("zero", "const_f64", TYPE_F64, attrs={"value": 0.0}),
            op("min0_hat", "min_f64", TYPE_F64, ("zero", "ell_hat")),
            op("stage1_pass", "lt_f64", TYPE_BOOL, ("log_u1", "min0_hat")),
        ),
        NativeTerminator("br_cond", ("stage1_pass",), {"true": "exact_candidate", "false": "reject_stage1"}),
    )
    reject1 = NativeBlock(
        "reject_stage1",
        (),
        NativeTerminator(
            "return", (), {"accepted": False, "reject_stage": 1, "exact_y_evaluated": False, "ell_R_valid": False}
        ),
    )
    exact_candidate = NativeBlock(
        "exact_candidate",
        (),
        NativeTerminator("br_cond", ("in.exact_x_valid",), {"true": "use_cache", "false": "compute_exact_x"}),
    )
    use_cache = NativeBlock(
        "use_cache",
        (op("cached_px", "load_input", TYPE_F64, attrs={"name": "exact_log_weight_x"}),),
        NativeTerminator("br", (), {"target": "after_exact_x"}),
    )
    compute_x = NativeBlock(
        "compute_exact_x",
        (
            op("exact_x", "call_graph", TYPE_C128, attrs={"graph": "exact", "binding": "x"}, region=REGION_EXACT),
            op(
                "computed_px",
                "log_positive_real",
                TYPE_F64,
                ("exact_x",),
                {"on_fail": "INVALID_EXACT_TARGET", "imag_atol": 1e-12},
                REGION_EXACT,
            ),
        ),
        NativeTerminator("br", (), {"target": "after_exact_x"}),
    )
    after = NativeBlock(
        "after_exact_x",
        (
            op(
                "log_px",
                "phi",
                TYPE_F64,
                (),
                {"preds": [["use_cache", "cached_px"], ["compute_exact_x", "computed_px"]]},
            ),
            op("exact_y", "call_graph", TYPE_C128, attrs={"graph": "exact", "binding": "y"}, region=REGION_EXACT),
            op(
                "log_py",
                "log_positive_real",
                TYPE_F64,
                ("exact_y",),
                {"on_fail": "INVALID_EXACT_TARGET", "imag_atol": 1e-12},
                REGION_EXACT,
            ),
            op("diff_pi", "sub_f64", TYPE_F64, ("log_py", "log_px")),
            op("ell_R", "add_f64", TYPE_F64, ("diff_pi", "log_q")),
            op("delta", "sub_f64", TYPE_F64, ("ell_R", "ell_hat")),
            op(
                "log_u2",
                "log_uniform",
                TYPE_F64,
                (),
                {"which": "u2", "on_fail": "INVALID_RANDOM_UNIFORM", "floor": 1e-300},
            ),
            op("min0_delta", "min_f64", TYPE_F64, ("zero", "delta")),
            op("stage2_pass", "lt_f64", TYPE_BOOL, ("log_u2", "min0_delta")),
        ),
        NativeTerminator("br_cond", ("stage2_pass",), {"true": "accept", "false": "reject_stage2"}),
    )
    reject2 = NativeBlock(
        "reject_stage2",
        (),
        NativeTerminator(
            "return", (), {"accepted": False, "reject_stage": 2, "exact_y_evaluated": True, "ell_R_valid": True}
        ),
    )
    accept = NativeBlock(
        "accept",
        (),
        NativeTerminator(
            "return", (), {"accepted": True, "reject_stage": 0, "exact_y_evaluated": True, "ell_R_valid": True}
        ),
    )
    return (entry, reject1, exact_candidate, use_cache, compute_x, after, reject2, accept)


def lower_native_kernel(kernel) -> NativeKernelIR:
    ir = kernel.ir
    layout = binding_layout_for(ir.family.model, ir.family.order)
    cheap_graph = lower_cheap_graph(kernel._cheap.dag)
    exact_dag = lower_diagramwise(ir, share=True)
    exact_graph = lower_exact_graph(exact_dag)
    source = {
        "cheap_policy_name": kernel.spec.cheap_policy_name,
        "design": kernel.spec.design,
        "diagram_ir_digest": kernel.spec.diagram_ir_digest,
        "diagram_ir_schema": kernel.spec.diagram_ir_schema,
        "exact_lowering": kernel.spec.exact_lowering,
        "numerical": dict(kernel.spec.numerical),
        "primitive_set": "native_kernel_v1",
        "proposal_kind": kernel.spec.proposal_kind,
        "target_policy_name": kernel.spec.target_policy_name,
    }
    return NativeKernelIR(
        SCHEMA_VERSION,
        IR_KIND,
        "diagram",
        source,
        layout,
        (cheap_graph, exact_graph),
        _control_blocks(source),
        "entry",
    )


def lower_native_score_kernel(*, proposal_kind: str = "provided_log_ratio") -> NativeKernelIR:
    """Design B CFG on explicit log-weights (finite-state revalidation)."""
    source = {
        "cheap_policy_name": "external_log_weights",
        "design": "B",
        "diagram_ir_digest": "none",
        "diagram_ir_schema": "1.0.0",
        "exact_lowering": "external_log_weights",
        "numerical": {"atol": 1e-12, "clip": "none", "rng_log_floor": 1e-300, "rtol": 1e-10},
        "proposal_kind": proposal_kind,
        "target_policy_name": "positive_real_F_v1",
    }

    def op(oid, name, type_, args=(), attrs=None):
        return NativeOp(oid, name, type_, tuple(args), dict(attrs or {}), REGION_CONTROL)

    layout = BindingLayout(SCHEMA_VERSION, "score", 0, ())
    entry = NativeBlock(
        "entry",
        (
            op("log_wx", "load_input", TYPE_F64, attrs={"name": "log_wx"}),
            op("log_wy", "load_input", TYPE_F64, attrs={"name": "log_wy"}),
            op("ell_hat", "sub_f64", TYPE_F64, ("log_wy", "log_wx")),
            op("log_q", "check_log_q", TYPE_F64, (), {"on_fail": "INVALID_PROPOSAL_RATIO", "kind": proposal_kind}),
            op(
                "log_u1",
                "log_uniform",
                TYPE_F64,
                (),
                {"which": "u1", "on_fail": "INVALID_RANDOM_UNIFORM", "floor": 1e-300},
            ),
            op("zero", "const_f64", TYPE_F64, attrs={"value": 0.0}),
            op("min0_hat", "min_f64", TYPE_F64, ("zero", "ell_hat")),
            op("stage1_pass", "lt_f64", TYPE_BOOL, ("log_u1", "min0_hat")),
        ),
        NativeTerminator("br_cond", ("stage1_pass",), {"true": "exact_candidate", "false": "reject_stage1"}),
    )
    reject1 = NativeBlock(
        "reject_stage1",
        (),
        NativeTerminator(
            "return", (), {"accepted": False, "reject_stage": 1, "exact_y_evaluated": False, "ell_R_valid": False}
        ),
    )
    exact_candidate = NativeBlock(
        "exact_candidate",
        (
            op("log_px", "load_input", TYPE_F64, attrs={"name": "log_px"}),
            op("log_py", "load_input", TYPE_F64, attrs={"name": "log_py"}),
            op("diff_pi", "sub_f64", TYPE_F64, ("log_py", "log_px")),
            op("ell_R", "add_f64", TYPE_F64, ("diff_pi", "log_q")),
            op("delta", "sub_f64", TYPE_F64, ("ell_R", "ell_hat")),
            op(
                "log_u2",
                "log_uniform",
                TYPE_F64,
                (),
                {"which": "u2", "on_fail": "INVALID_RANDOM_UNIFORM", "floor": 1e-300},
            ),
            op("min0_delta", "min_f64", TYPE_F64, ("zero", "delta")),
            op("stage2_pass", "lt_f64", TYPE_BOOL, ("log_u2", "min0_delta")),
        ),
        NativeTerminator("br_cond", ("stage2_pass",), {"true": "accept", "false": "reject_stage2"}),
    )
    reject2 = NativeBlock(
        "reject_stage2",
        (),
        NativeTerminator(
            "return", (), {"accepted": False, "reject_stage": 2, "exact_y_evaluated": True, "ell_R_valid": True}
        ),
    )
    accept = NativeBlock(
        "accept",
        (),
        NativeTerminator(
            "return", (), {"accepted": True, "reject_stage": 0, "exact_y_evaluated": True, "ell_R_valid": True}
        ),
    )
    return NativeKernelIR(
        SCHEMA_VERSION,
        IR_KIND,
        "score",
        source,
        layout,
        (),
        (entry, reject1, exact_candidate, reject2, accept),
        "entry",
    )
