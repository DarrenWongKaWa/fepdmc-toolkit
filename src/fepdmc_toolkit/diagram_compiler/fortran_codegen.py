"""NativeKernelIR -> deterministic Fortran source. Consumes NativeKernelIR only."""

from __future__ import annotations

import hashlib
import json
import re

from .native_ir import (
    INVALID_IR,
    TYPE_BOOL,
    TYPE_C128,
    TYPE_C128_M2,
    TYPE_F64,
    TYPE_I64,
    NativeIRError,
    NativeKernelIR,
)

BACKEND = "fortran_v1"
RUNTIME = "native_kernel_v1"
EPILOGUE_LABEL = 900

_KNOWN_LABELS = {
    "entry": 100,
    "reject_stage1": 200,
    "exact_candidate": 300,
    "use_cache": 310,
    "compute_exact_x": 320,
    "after_exact_x": 330,
    "reject_stage2": 400,
    "accept": 500,
}

SUPPORTED_OPS = frozenset(
    {
        "const_f64",
        "const_c128",
        "const_c128_m2",
        "load",
        "load_index",
        "load_input",
        "add_f64",
        "sub_f64",
        "mul_f64",
        "min_f64",
        "lt_f64",
        "add_c128",
        "mul_c128",
        "add_m2",
        "scale_m2",
        "matmul",
        "trace",
        "prim.electron_scalar",
        "prim.electron_twoband",
        "prim.phonon",
        "prim.prefactor",
        "prim.vertex_sigmaz",
        "log_positive_real",
        "log_uniform",
        "check_log_q",
        "call_graph",
        "phi",
    }
)

SUPPORTED_TERMS = frozenset({"br", "br_cond", "return"})

_INPUT_DUMMY = {
    "exact_log_weight_x": "exact_log_weight_x",
    "log_wx": "log_wx_in",
    "log_wy": "log_wy_in",
    "log_px": "log_px_in",
    "log_py": "log_py_in",
    "log_q": "log_q_reverse_minus_forward",
}

_OUTPUT_DUMMIES = (
    "accepted",
    "reject_stage",
    "ell_hat",
    "ell_R_valid",
    "ell_R",
    "exact_y_evaluated",
    "status",
)


def native_ir_digest(nir: NativeKernelIR) -> str:
    payload = json.dumps(nir.to_dict(), sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(payload).hexdigest()


def _ident(name: str) -> str:
    text = name.replace(".", "_").replace("-", "_")
    if not text or (not text[0].isalpha() and text[0] != "_"):
        text = "v_" + text
    if not re.fullmatch(r"[A-Za-z_][A-Za-z0-9_]*", text):
        raise NativeIRError(INVALID_IR, f"cannot emit identifier for {name}")
    return text


def _f64_lit(value: float) -> str:
    return format(float(value), ".17e") + "_real64"


def _i64_lit(value: int) -> str:
    return f"{int(value)}_int64"


def _labels_for(nir: NativeKernelIR) -> dict[str, int]:
    used = set(_KNOWN_LABELS.values())
    used.add(EPILOGUE_LABEL)
    labels: dict[str, int] = {}
    for block in nir.blocks:
        if block.id in _KNOWN_LABELS:
            labels[block.id] = _KNOWN_LABELS[block.id]
            continue
        digest = hashlib.sha256(block.id.encode("utf-8")).hexdigest()
        lab = 600 + (int(digest[:4], 16) % 250)
        while lab in used:
            lab += 1
            if lab > 899:
                lab = 600
        used.add(lab)
        labels[block.id] = lab
    return labels


def _decl(type_name: str, name: str) -> str:
    ident = _ident(name)
    if type_name == TYPE_F64:
        return f"real(real64) :: {ident}"
    if type_name == TYPE_C128:
        return f"complex(real64) :: {ident}"
    if type_name == TYPE_C128_M2:
        return f"complex(real64) :: {ident}(2, 2)"
    if type_name == TYPE_BOOL:
        return f"logical :: {ident}"
    if type_name == TYPE_I64:
        return f"integer(int64) :: {ident}"
    raise NativeIRError(INVALID_IR, f"unsupported type {type_name}")


def _layout_args(layout, prefix: str) -> list[tuple[str, str, tuple]]:
    return [(f"{prefix}_{field.name}", field.type, field.shape) for field in layout.fields]


def _wrap_args(names: list[str], indent: str = "    ") -> str:
    if not names:
        return ""
    lines: list[str] = []
    current = indent
    for i, name in enumerate(names):
        piece = name + (", " if i < len(names) - 1 else "")
        if current != indent and len(current) + len(piece) > 88:
            lines.append(current + "&")
            current = indent + "  " + piece
        else:
            current += piece
    lines.append(current)
    return "\n".join(lines)


def _kernel_arg_names(nir: NativeKernelIR) -> list[str]:
    names: list[str] = []
    if nir.kind == "score":
        names.extend(["log_wx_in", "log_wy_in", "log_px_in", "log_py_in"])
    else:
        for prefix in ("x", "y"):
            names.extend(_ident(name) for name, _t, _s in _layout_args(nir.layout, prefix))
        names.extend(["log_q_reverse_minus_forward", "u1", "u2", "exact_x_valid", "exact_log_weight_x"])
        return names + list(_OUTPUT_DUMMIES)
    names.extend(["log_q_reverse_minus_forward", "u1", "u2"])
    return names + list(_OUTPUT_DUMMIES)


def _dummy_decls(nir: NativeKernelIR) -> list[str]:
    lines: list[str] = []
    if nir.kind == "score":
        lines.append("real(real64), intent(in) :: log_wx_in, log_wy_in, log_px_in, log_py_in")
        lines.append("real(real64), intent(in) :: log_q_reverse_minus_forward, u1, u2")
    else:
        for prefix in ("x", "y"):
            for name, type_name, shape in _layout_args(nir.layout, prefix):
                ident = _ident(name)
                if shape == ():
                    if type_name == TYPE_F64:
                        lines.append(f"real(real64), intent(in) :: {ident}")
                    elif type_name == TYPE_I64:
                        lines.append(f"integer(int64), intent(in) :: {ident}")
                    else:
                        raise NativeIRError(INVALID_IR, f"unsupported dummy type {type_name}")
                else:
                    extents = ", ".join(str(int(n)) for n in shape)
                    lines.append(f"real(real64), intent(in) :: {ident}({extents})")
        lines.append("real(real64), intent(in) :: log_q_reverse_minus_forward, u1, u2")
        lines.append("logical, intent(in) :: exact_x_valid")
        lines.append("real(real64), intent(in) :: exact_log_weight_x")
    lines.append("logical, intent(out) :: accepted, ell_R_valid, exact_y_evaluated")
    lines.append("integer(int64), intent(out) :: reject_stage, status")
    lines.append("real(real64), intent(out) :: ell_hat, ell_R")
    return lines


def _graph_arg_names(layout) -> list[str]:
    return [_ident(name) for name, _t, _s in _layout_args(layout, "b")] + ["graph_value"]


def _graph_dummy_decls(layout) -> list[str]:
    lines: list[str] = []
    for name, type_name, shape in _layout_args(layout, "b"):
        ident = _ident(name)
        if shape == ():
            if type_name == TYPE_F64:
                lines.append(f"real(real64), intent(in) :: {ident}")
            elif type_name == TYPE_I64:
                lines.append(f"integer(int64), intent(in) :: {ident}")
            else:
                raise NativeIRError(INVALID_IR, f"unsupported graph dummy {type_name}")
        else:
            extents = ", ".join(str(int(n)) for n in shape)
            lines.append(f"real(real64), intent(in) :: {ident}({extents})")
    lines.append("complex(real64), intent(out) :: graph_value")
    return lines


def _graph_call_args(layout, prefix: str) -> list[str]:
    return [_ident(name) for name, _t, _s in _layout_args(layout, prefix)]


def _emit_op(op, *, graph_binding: str | None, layout, allow_status_goto: bool) -> list[str]:
    dest = _ident(op.id)
    args = [_ident(a) for a in op.args]
    name = op.op
    if name not in SUPPORTED_OPS:
        raise NativeIRError(INVALID_IR, f"unsupported opcode {name}")
    if name == "const_f64":
        return [f"{dest} = {_f64_lit(op.attrs['value'])}"]
    if name == "const_c128":
        return [f"{dest} = cmplx({_f64_lit(op.attrs['re'])}, {_f64_lit(op.attrs['im'])}, kind=real64)"]
    if name == "const_c128_m2":
        lines = []
        for i, row in enumerate(op.attrs["value"]):
            for j, pair in enumerate(row):
                re, im = pair
                lines.append(f"{dest}({i + 1}, {j + 1}) = cmplx({_f64_lit(re)}, {_f64_lit(im)}, kind=real64)")
        return lines
    if name == "load":
        field = op.attrs["field"]
        prefix = "b" if graph_binding else op.attrs.get("binding", "x")
        src = _ident(f"{prefix}_{field}")
        if field == "L":
            return [f"{dest} = real({src}, real64)"]
        return [f"{dest} = {src}"]
    if name == "load_index":
        field = op.attrs["field"]
        index = int(op.attrs["index"]) + 1
        prefix = "b" if graph_binding else op.attrs.get("binding", "x")
        src = _ident(f"{prefix}_{field}")
        return [f"{dest} = {src}({index})"]
    if name == "load_input":
        key = op.attrs["name"]
        if key not in _INPUT_DUMMY:
            raise NativeIRError(INVALID_IR, f"unsupported load_input {key}")
        return [f"{dest} = {_INPUT_DUMMY[key]}"]
    if name == "add_f64":
        return [f"{dest} = {args[0]} + {args[1]}"]
    if name == "sub_f64":
        return [f"{dest} = {args[0]} - {args[1]}"]
    if name == "mul_f64":
        return [f"{dest} = {args[0]} * {args[1]}"]
    if name == "min_f64":
        return [f"{dest} = min({args[0]}, {args[1]})"]
    if name == "lt_f64":
        return [f"{dest} = ({args[0]} < {args[1]})"]
    if name == "add_c128":
        return [f"{dest} = {args[0]} + {args[1]}"]
    if name == "mul_c128":
        return [f"{dest} = {args[0]} * {args[1]}"]
    if name == "add_m2":
        return [f"call nk_add_m2({args[0]}, {args[1]}, {dest})"]
    if name == "scale_m2":
        return [f"call nk_scale_m2({args[0]}, {args[1]}, {dest})"]
    if name == "matmul":
        return [f"call nk_matmul({args[0]}, {args[1]}, {dest})"]
    if name == "trace":
        return [f"call nk_trace_counted({args[0]}, {dest})"]
    if name == "prim.electron_scalar":
        return [f"{dest} = nk_electron_scalar({args[0]}, {args[1]}, {args[2]})"]
    if name == "prim.electron_twoband":
        return [f"call nk_electron_twoband({args[0]}, {args[1]}, {args[2]}, {args[3]}, {args[4]}, {dest})"]
    if name == "prim.phonon":
        return [f"{dest} = nk_phonon({args[0]}, {args[1]})"]
    if name == "prim.prefactor":
        n = int(op.attrs["n"])
        return [f"{dest} = nk_prefactor({args[0]}, {args[1]}, {_i64_lit(n)})"]
    if name == "prim.vertex_sigmaz":
        return [f"call nk_vertex_sigmaz({args[0]}, {args[1]}, {dest})"]
    if name == "log_positive_real":
        if not allow_status_goto:
            raise NativeIRError(INVALID_IR, "log_positive_real is control-only")
        fail = op.attrs.get("on_fail", "INVALID_EXACT_TARGET")
        code = {
            "INVALID_CHEAP_WEIGHT": "NK_INVALID_CHEAP_WEIGHT",
            "INVALID_EXACT_TARGET": "NK_INVALID_EXACT_TARGET",
        }.get(fail)
        if code is None:
            raise NativeIRError(INVALID_IR, f"unknown on_fail {fail}")
        atol = float(op.attrs.get("imag_atol", 1e-12))
        return [
            f"call nk_log_positive_real({args[0]}, {_f64_lit(atol)}, {dest}, status, {code})",
            f"if (status /= NK_OK) goto {EPILOGUE_LABEL}",
        ]
    if name == "log_uniform":
        if not allow_status_goto:
            raise NativeIRError(INVALID_IR, "log_uniform is control-only")
        which = op.attrs["which"]
        if which not in ("u1", "u2"):
            raise NativeIRError(INVALID_IR, f"unsupported uniform {which}")
        floor = float(op.attrs.get("floor", 1e-300))
        return [
            f"call nk_log_uniform({which}, {_f64_lit(floor)}, {dest}, status)",
            f"if (status /= NK_OK) goto {EPILOGUE_LABEL}",
        ]
    if name == "check_log_q":
        if not allow_status_goto:
            raise NativeIRError(INVALID_IR, "check_log_q is control-only")
        kind = op.attrs.get("kind", "provided_log_ratio")
        flag = "1_int64" if kind == "symmetric" else "0_int64"
        return [
            f"call nk_check_log_q(log_q_reverse_minus_forward, {flag}, {dest}, status)",
            f"if (status /= NK_OK) goto {EPILOGUE_LABEL}",
        ]
    if name == "call_graph":
        graph = op.attrs["graph"]
        binding = op.attrs["binding"]
        call_args = _graph_call_args(layout, binding) + [dest]
        joined = ", ".join(call_args)
        return [f"call eval_graph_{_ident(graph)}({joined})"]
    if name == "phi":
        return []
    raise NativeIRError(INVALID_IR, f"unsupported opcode {name}")


def _collect_kernel_decls(nir: NativeKernelIR) -> list[str]:
    seen: dict[str, str] = {}
    for block in nir.blocks:
        for op in block.ops:
            seen[op.id] = op.type
    skip = set(_OUTPUT_DUMMIES)
    lines = []
    for name in sorted(seen):
        if name in skip:
            continue
        lines.append(_decl(seen[name], name))
    return lines


def _emit_graph(graph, layout) -> str:
    args = _graph_arg_names(layout)
    header = f"  subroutine eval_graph_{_ident(graph.id)}(&\n" + _wrap_args(args, "      ") + ")"
    lines = [header]
    for dummy in _graph_dummy_decls(layout):
        lines.append("    " + dummy)
    seen = {op.id: op.type for op in graph.ops}
    for name in sorted(seen):
        lines.append("    " + _decl(seen[name], name))
    if graph.region == "exact":
        lines.append("    call nk_bump_exact_graph()")
    elif graph.region == "cheap":
        lines.append("    call nk_bump_cheap_graph()")
    for op in graph.ops:
        for stmt in _emit_op(op, graph_binding="b", layout=layout, allow_status_goto=False):
            lines.append("    " + stmt)
    lines.append(f"    graph_value = {_ident(graph.output)}")
    lines.append(f"  end subroutine eval_graph_{_ident(graph.id)}")
    return "\n".join(lines)


def _phi_copies(nir: NativeKernelIR, src_block: str, dst_block: str) -> list[str]:
    lines = []
    target = nir.block(dst_block)
    for op in target.ops:
        if op.op != "phi":
            continue
        for pred, value in op.attrs["preds"]:
            if pred == src_block:
                lines.append(f"{_ident(op.id)} = {_ident(value)}")
    return lines


def _emit_goto(nir, labels, src: str, dst: str) -> list[str]:
    lines = _phi_copies(nir, src, dst)
    lines.append(f"goto {labels[dst]}")
    return lines


def _logical_lit(flag: bool) -> str:
    return ".true." if flag else ".false."


def emit_fortran(nir: NativeKernelIR, backend: str = BACKEND) -> str:
    if backend != BACKEND:
        raise NativeIRError(INVALID_IR, f"unsupported backend {backend}")
    if nir.ir_kind != "NativeKernelIR":
        raise NativeIRError(INVALID_IR, "backend consumes NativeKernelIR only")
    prim = str(nir.source.get("primitive_set", RUNTIME))
    if prim != RUNTIME:
        raise NativeIRError(INVALID_IR, f"unsupported primitive_set {prim}")
    for graph in nir.graphs:
        for op in graph.ops:
            if op.op not in SUPPORTED_OPS:
                raise NativeIRError(INVALID_IR, f"unsupported opcode {op.op}")
    for block in nir.blocks:
        for op in block.ops:
            if op.op not in SUPPORTED_OPS:
                raise NativeIRError(INVALID_IR, f"unsupported opcode {op.op}")
        if block.term.op not in SUPPORTED_TERMS:
            raise NativeIRError(INVALID_IR, f"unsupported terminator {block.term.op}")
    digest = native_ir_digest(nir)
    labels = _labels_for(nir)
    layout = nir.layout
    arg_names = _kernel_arg_names(nir)
    lines = [
        f"! backend={BACKEND} primitive_set={prim}",
        f"! schema_version={nir.schema_version} native_ir_digest={digest}",
        f"! ir_kind={nir.kind} design={nir.source.get('design') or ''} "
        f"cheap_policy={nir.source.get('cheap_policy_name') or ''}",
        "! cache_contract=caller_guarantees_exact_x_valid_belongs_to_x",
        "module generated_da_kernel",
        "  use, intrinsic :: iso_fortran_env, only: real64, int64",
        "  use native_kernel_runtime_v1",
        "  implicit none",
        f"  character(len=*), parameter :: NK_BACKEND = '{BACKEND}'",
        f"  character(len=*), parameter :: NK_PRIMITIVE_SET = '{prim}'",
        f"  character(len=*), parameter :: NK_IR_DIGEST = '{digest}'",
        f"  integer, parameter :: NK_SCHEMA_VERSION = {int(nir.schema_version)}",
        "contains",
    ]
    for graph in nir.graphs:
        lines.append(_emit_graph(graph, layout))
    lines.append("  subroutine evaluate_da_kernel(&")
    lines.append(_wrap_args(arg_names, "      ") + ")")
    for dummy in _dummy_decls(nir):
        lines.append("    " + dummy)
    for decl in _collect_kernel_decls(nir):
        lines.append("    " + decl)
    lines.append("    call nk_reset_counters()")
    lines.append("    status = NK_OK")
    lines.append("    accepted = .false.")
    lines.append("    reject_stage = 0_int64")
    lines.append("    ell_R_valid = .false.")
    lines.append("    exact_y_evaluated = .false.")
    lines.append("    ell_hat = 0.0_real64")
    lines.append("    ell_R = 0.0_real64")
    lines.append(f"    goto {labels[nir.entry]}")
    for block in nir.blocks:
        lines.append(f"    ! block {block.id}")
        lines.append(f"{labels[block.id]} continue")
        for op in block.ops:
            if op.op == "phi":
                continue
            for stmt in _emit_op(op, graph_binding=None, layout=layout, allow_status_goto=True):
                lines.append("    " + stmt)
        term = block.term
        if term.op == "br":
            for stmt in _emit_goto(nir, labels, block.id, str(term.attrs["target"])):
                lines.append("    " + stmt)
        elif term.op == "br_cond":
            cond = term.args[0]
            if cond == "in.exact_x_valid":
                cond_expr = "exact_x_valid"
            else:
                cond_expr = _ident(cond)
            true_b = str(term.attrs["true"])
            false_b = str(term.attrs["false"])
            lines.append(f"    if ({cond_expr}) then")
            for stmt in _emit_goto(nir, labels, block.id, true_b):
                lines.append("      " + stmt)
            lines.append("    else")
            for stmt in _emit_goto(nir, labels, block.id, false_b):
                lines.append("      " + stmt)
            lines.append("    end if")
        elif term.op == "return":
            lines.append(f"    accepted = {_logical_lit(bool(term.attrs['accepted']))}")
            lines.append(f"    reject_stage = {_i64_lit(int(term.attrs['reject_stage']))}")
            lines.append(f"    exact_y_evaluated = {_logical_lit(bool(term.attrs['exact_y_evaluated']))}")
            lines.append(f"    ell_R_valid = {_logical_lit(bool(term.attrs['ell_R_valid']))}")
            lines.append(f"    goto {EPILOGUE_LABEL}")
        else:
            raise NativeIRError(INVALID_IR, f"unsupported terminator {term.op}")
    lines.append(f"{EPILOGUE_LABEL} continue")
    lines.append("    if (status /= NK_OK) then")
    lines.append("      accepted = .false.")
    lines.append("      reject_stage = 0_int64")
    lines.append("      ell_R_valid = .false.")
    lines.append("      exact_y_evaluated = .false.")
    lines.append("    end if")
    lines.append("  end subroutine evaluate_da_kernel")
    lines.append("end module generated_da_kernel")
    source = "\n".join(lines) + "\n"
    validate_generated_fortran(source, nir)
    return source


def emit_driver(nir: NativeKernelIR) -> str:
    arg_names = _kernel_arg_names(nir)
    lines = [
        "program nk_eval_driver",
        "  use, intrinsic :: iso_fortran_env, only: real64, int64",
        "  use native_kernel_runtime_v1",
        "  use generated_da_kernel",
        "  implicit none",
        "  logical :: accepted, ell_R_valid, exact_y_evaluated",
        "  integer(int64) :: reject_stage, status, accepted_i, ell_R_valid_i, exact_y_i",
        "  real(real64) :: ell_hat, ell_R",
        "  real(real64) :: log_q_reverse_minus_forward, u1, u2",
    ]
    if nir.kind == "score":
        lines.append("  real(real64) :: log_wx_in, log_wy_in, log_px_in, log_py_in")
        lines.append("  read(*,*) log_wx_in, log_wy_in, log_px_in, log_py_in")
        lines.append("  read(*,*) log_q_reverse_minus_forward, u1, u2")
    else:
        lines.append("  integer(int64) :: exact_x_valid_i")
        lines.append("  logical :: exact_x_valid")
        lines.append("  real(real64) :: exact_log_weight_x")
        for prefix in ("x", "y"):
            for name, type_name, shape in _layout_args(nir.layout, prefix):
                ident = _ident(name)
                if shape == ():
                    if type_name == TYPE_F64:
                        lines.append(f"  real(real64) :: {ident}")
                    else:
                        lines.append(f"  integer(int64) :: {ident}")
                else:
                    extents = ", ".join(str(int(n)) for n in shape)
                    lines.append(f"  real(real64) :: {ident}({extents})")
        for prefix in ("x", "y"):
            for name, _type, _shape in _layout_args(nir.layout, prefix):
                lines.append(f"  read(*,*) {_ident(name)}")
        lines.append("  read(*,*) log_q_reverse_minus_forward, u1, u2")
        lines.append("  read(*,*) exact_x_valid_i, exact_log_weight_x")
        lines.append("  exact_x_valid = (exact_x_valid_i /= 0_int64)")
    lines.append("  call evaluate_da_kernel(&")
    lines.append(_wrap_args(arg_names, "    ") + ")")
    lines.append("  accepted_i = merge(1_int64, 0_int64, accepted)")
    lines.append("  ell_R_valid_i = merge(1_int64, 0_int64, ell_R_valid)")
    lines.append("  exact_y_i = merge(1_int64, 0_int64, exact_y_evaluated)")
    lines.append("  write(*,'(I0,1X,I0,1X,I0,1X,ES24.16E3,1X,I0,1X,ES24.16E3,1X,I0,1X,I0,1X,I0,1X,I0,1X,I0,1X,I0)') &")
    lines.append("    status, accepted_i, reject_stage, ell_hat, ell_R_valid_i, ell_R, exact_y_i, &")
    lines.append("    nk_n_electron_twoband, nk_n_vertex_sigmaz, nk_n_matmul, nk_n_trace, nk_n_exact_graph")
    lines.append("end program nk_eval_driver")
    return "\n".join(lines) + "\n"


_EXACT_PRIMITIVE_MARKERS = (
    "eval_graph_exact",
    "nk_electron_twoband",
    "nk_vertex_sigmaz",
    "nk_matmul",
    "nk_trace_counted",
)


def validate_generated_fortran(source: str, nir: NativeKernelIR) -> None:
    if "\r" in source:
        raise NativeIRError(INVALID_IR, "generated Fortran contains CR")
    if "None" in source or "True" in source or "False" in source:
        raise NativeIRError(INVALID_IR, "Python literal leaked into Fortran")
    if "TODO" in source or "placeholder" in source.lower():
        raise NativeIRError(INVALID_IR, "unresolved placeholder in Fortran")
    if BACKEND not in source or RUNTIME not in source:
        raise NativeIRError(INVALID_IR, "missing backend/runtime identity")
    digest = native_ir_digest(nir)
    if digest not in source:
        raise NativeIRError(INVALID_IR, "missing NativeKernelIR digest")
    labels = _labels_for(nir)
    for block in nir.blocks:
        token = f"{labels[block.id]} continue"
        if token not in source:
            raise NativeIRError(INVALID_IR, f"missing block label {block.id}")
        for pred, _value in (item for op in block.ops if op.op == "phi" for item in op.attrs["preds"]):
            if pred not in labels:
                raise NativeIRError(INVALID_IR, f"phi pred {pred} missing")
    kernel = source.split("subroutine evaluate_da_kernel", 1)[1]
    kernel_body = kernel.split("end subroutine evaluate_da_kernel", 1)[0]
    if nir.kind == "diagram":
        entry_lab = labels["entry"]
        reject1 = labels.get("reject_stage1")
        exact_lab = labels.get("exact_candidate")
        if exact_lab is not None and reject1 is not None:
            entry_chunk = kernel_body.split(f"{entry_lab} continue", 1)[1]
            entry_chunk = entry_chunk.split(f"{reject1} continue", 1)[0]
            for marker in _EXACT_PRIMITIVE_MARKERS:
                if marker in entry_chunk:
                    raise NativeIRError(INVALID_IR, f"exact work in Stage-1 chunk: {marker}")
            if f"{reject1} continue" in kernel_body:
                rchunk = kernel_body.split(f"{reject1} continue", 1)[1]
                rchunk = rchunk.split(f"{exact_lab} continue", 1)[0]
                for marker in _EXACT_PRIMITIVE_MARKERS:
                    if marker in rchunk:
                        raise NativeIRError(INVALID_IR, f"exact work on Stage-1 reject: {marker}")
    kernel_text = kernel_body
    for block in nir.blocks:
        for op in block.ops:
            ident = _ident(op.id)
            if op.op == "phi":
                copies = 0
                for _pred, value in op.attrs["preds"]:
                    needle = f"{ident} = {_ident(value)}"
                    copies += kernel_text.count(needle)
                if copies != len(op.attrs["preds"]):
                    raise NativeIRError(INVALID_IR, f"phi {op.id} not copied on every edge")
                continue
            if ident not in kernel_text:
                raise NativeIRError(INVALID_IR, f"SSA {op.id} not emitted")
    for graph in nir.graphs:
        sub = f"subroutine eval_graph_{_ident(graph.id)}"
        if sub not in source:
            raise NativeIRError(INVALID_IR, f"missing graph {graph.id}")
        for op in graph.ops:
            ident = _ident(op.id)
            if ident not in source:
                raise NativeIRError(INVALID_IR, f"graph SSA {op.id} missing")
    if "error stop" in kernel_body:
        raise NativeIRError(INVALID_IR, "generated kernel must not error stop")
    for field in nir.layout.fields:
        for prefix in ("x", "y"):
            ident = _ident(f"{prefix}_{field.name}")
            if ident not in source and nir.kind == "diagram":
                raise NativeIRError(INVALID_IR, f"missing ABI field {ident}")
    if "NK_INVALID_CHEAP_WEIGHT" not in source and nir.kind == "diagram":
        raise NativeIRError(INVALID_IR, "cheap target validation missing")
    if "nk_check_log_q" not in source:
        raise NativeIRError(INVALID_IR, "proposal ratio check missing")
