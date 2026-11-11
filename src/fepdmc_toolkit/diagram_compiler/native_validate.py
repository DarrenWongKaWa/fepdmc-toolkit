"""Static validator and exact-laziness reachability for NativeKernelIR."""

from __future__ import annotations

from .native_ir import INVALID_IR, REGION_EXACT, NativeIRError, NativeKernelIR

_EXACT_MARKERS = frozenset(
    {"prim.electron_twoband", "prim.vertex_sigmaz", "matmul", "trace", "add_m2", "scale_m2", "const_c128_m2"}
)


def _successors(block) -> tuple[str, ...]:
    term = block.term
    if term.op == "br":
        return (str(term.attrs["target"]),)
    if term.op == "br_cond":
        return (str(term.attrs["true"]), str(term.attrs["false"]))
    if term.op == "return":
        return ()
    raise NativeIRError(INVALID_IR, f"unknown terminator {term.op}")


def _defined_ids(nir: NativeKernelIR) -> set[str]:
    ids = {"in.exact_x_valid", "in.exact_log_weight_x", "in.log_q", "in.u1", "in.u2"}
    for graph in nir.graphs:
        seen: set[str] = set()
        for op in graph.ops:
            if op.id in seen or op.id in ids:
                raise NativeIRError(INVALID_IR, f"duplicate id {op.id}")
            seen.add(op.id)
        if graph.output not in seen:
            raise NativeIRError(INVALID_IR, f"graph {graph.id} output undefined")
    for block in nir.blocks:
        for op in block.ops:
            if op.id in ids:
                raise NativeIRError(INVALID_IR, f"duplicate id {op.id}")
            ids.add(op.id)
    return ids


def cfg_successors(nir: NativeKernelIR) -> dict[str, tuple[str, ...]]:
    return {block.id: _successors(block) for block in nir.blocks}


def reachable(nir: NativeKernelIR, start: str) -> set[str]:
    succ = cfg_successors(nir)
    seen = set()
    stack = [start]
    while stack:
        name = stack.pop()
        if name in seen:
            continue
        seen.add(name)
        stack.extend(succ.get(name, ()))
    return seen


def stage1_reject_blocks(nir: NativeKernelIR) -> set[str]:
    """Blocks reachable from entry without taking the Stage-1-pass edge."""
    entry = nir.block(nir.entry)
    if entry.term.op != "br_cond":
        raise NativeIRError(INVALID_IR, "entry must branch on Stage 1")
    reject = str(entry.term.attrs["false"])
    return {nir.entry, reject} | reachable(nir, reject)


def exact_candidate_ops_on_reject_path(nir: NativeKernelIR) -> list[str]:
    banned = []
    reject_set = stage1_reject_blocks(nir)
    for block in nir.blocks:
        if block.id not in reject_set:
            continue
        for op in block.ops:
            if op.region == REGION_EXACT or op.op in _EXACT_MARKERS:
                banned.append(f"{block.id}:{op.id}")
            if op.op == "call_graph" and op.attrs.get("graph") == "exact":
                banned.append(f"{block.id}:{op.id}")
    return banned


def validate_native_ir(nir: NativeKernelIR) -> None:
    if nir.ir_kind != "NativeKernelIR" or nir.schema_version != 1:
        raise NativeIRError(INVALID_IR, "bad NativeKernelIR header")
    block_ids = [block.id for block in nir.blocks]
    if len(block_ids) != len(set(block_ids)):
        raise NativeIRError(INVALID_IR, "duplicate block id")
    if nir.entry not in block_ids:
        raise NativeIRError(INVALID_IR, "missing entry")
    _defined_ids(nir)
    succ = cfg_successors(nir)
    for name, targets in succ.items():
        for target in targets:
            if target not in block_ids:
                raise NativeIRError(INVALID_IR, f"invalid branch target {target} from {name}")
    for block in nir.blocks:
        for op in block.ops:
            if op.op == "phi":
                continue
            for arg in op.args:
                if arg.startswith("in."):
                    continue
    banned = exact_candidate_ops_on_reject_path(nir)
    if banned:
        raise NativeIRError(INVALID_IR, f"exact ops on Stage-1 reject path: {banned}")
    if nir.kind == "diagram":
        names = {graph.id for graph in nir.graphs}
        if "cheap" not in names or "exact" not in names:
            raise NativeIRError(INVALID_IR, "diagram kernel needs cheap and exact graphs")
        cheap = nir.graph("cheap")
        for op in cheap.ops:
            if op.op in _EXACT_MARKERS or op.region == REGION_EXACT:
                raise NativeIRError(INVALID_IR, "exact primitive in cheap graph")
