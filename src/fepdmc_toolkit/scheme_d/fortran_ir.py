"""Actual Graph -> Fortran -> shared-library execution with runtime bindings.

Straight-line emitted IR, including noncommuting matrices and scale annotations.
This standalone backend is NOT a perturbo overlay. One instance owns reusable
scratch/output buffers; use separate instances for simultaneous threads.
"""

from __future__ import annotations

import ctypes
import hashlib
import shutil
import subprocess
import time
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from .exec import _as_matrix
from .ir import Graph, Layout, Op


def _number(z):
    z = complex(z)
    if not np.isfinite(z):
        raise ValueError("nonfinite constant/annotation")

    def real(v):
        return format(float(v), ".17e").replace("e", "d")

    return f"cmplx({real(z.real)},{real(z.imag)},kind=c_double)"


@dataclass
class Emitted:
    source: str
    variables: dict
    outputs: dict
    n_input: int
    n_work: int
    n_output: int
    live_nodes: int


def emit_fortran(graph: Graph) -> Emitted:
    live = set()
    if not graph.outputs:
        raise ValueError("no graph outputs")
    pending = list(graph.outputs.values())
    while pending:
        i = pending.pop()
        if not 0 <= i < len(graph.nodes):
            raise ValueError("invalid output/dependency")
        if i in live:
            continue
        live.add(i)
        for j in graph.nodes[i].inputs:
            if j >= i:
                raise ValueError("graph is not topologically ordered")
            pending.append(j)
    variables = {}
    offsets = {}
    n_input = 0
    n_work = 0
    for i, node in enumerate(graph.nodes):
        if i not in live:
            continue
        if len(node.shape) != 2 or any(not isinstance(k, int) or k <= 0 for k in node.shape):
            raise ValueError("backend requires positive rank-two tensor shapes")
        size = int(np.prod(node.shape))
        offsets[i] = (n_work + 1, n_work + size)
        n_work += size
        if node.op == Op.INPUT.value:
            spec = (node.shape, node.layout, node.binding, node.material_id, node.gauge_id)
            if node.var_name in variables:
                if variables[node.var_name]["spec"] != spec:
                    raise ValueError("ambiguous input binding")
            else:
                variables[node.var_name] = {"offset": n_input, "size": size, "spec": spec}
                n_input += size
    outputs = {}
    n_output = 0
    for name, i in graph.outputs.items():
        size = int(np.prod(graph.nodes[i].shape))
        outputs[name] = {"offset": n_output, "size": size, "shape": graph.nodes[i].shape}
        n_output += size
    # Apply identical liveness allocation to flat, generic-CSE and typed graphs.
    # Dead intermediates are not kept merely to handicap the flat baseline.
    last = {i: i for i in live}
    for i in sorted(live):
        for j in graph.nodes[i].inputs:
            last[j] = max(last[j], i)
    for i in graph.outputs.values():
        last[i] = len(graph.nodes)
    expires = {}
    for i, end in last.items():
        expires.setdefault(end, []).append(i)
    free = []
    offsets = {}
    n_work = 0
    for i in sorted(live):
        size = int(np.prod(graph.nodes[i].shape))
        fits = [(length, j) for j, (_start, length) in enumerate(free) if length >= size]
        if fits:
            _, j = min(fits)
            start, length = free.pop(j)
            if length > size:
                free.append((start + size, length - size))
        else:
            start = n_work + 1
            n_work += size
        offsets[i] = (start, start + size - 1)
        for old in expires.get(i, []):
            if old != i:
                a, z = offsets[old]
                free.append((a, z - a + 1))

    def flat(i):
        a, b = offsets[i]
        return f"w({a}:{b})"

    def matrix(i):
        r, c = graph.nodes[i].shape
        return f"reshape({flat(i)},[{r},{c}])"

    code = [
        "module scheme_d_generated",
        "use iso_c_binding",
        "use, intrinsic :: ieee_arithmetic",
        "implicit none",
        "contains",
        'subroutine evaluate(x,w,y,status) bind(C,name="scheme_d_eval")',
        f"complex(c_double_complex), intent(in) :: x({max(1, n_input)})",
        f"complex(c_double_complex), intent(inout) :: w({max(1, n_work)})",
        f"complex(c_double_complex), intent(out) :: y({n_output})",
        "integer(c_int), intent(out) :: status",
        "integer :: ir,ic,it",
        "real(c_double) :: nm",
        "status=0",
    ]
    for i, node in enumerate(graph.nodes):
        if i not in live:
            continue
        dest = flat(i)
        r, c = node.shape
        size = r * c
        op = node.op
        ins = node.inputs
        if op == Op.INPUT.value:
            k = variables[node.var_name]["offset"]
            expr = f"x({k + 1}:{k + size})"
        elif op == Op.CONST.value:
            expr = _number(complex(dict(node.attrs)["value"]))
        elif op in {Op.ADD.value, Op.MUL.value}:
            symbol = "+" if op == Op.ADD.value else "*"
            expr = f"{flat(ins[0])} {symbol} {flat(ins[1])}"
        elif op == Op.MATMUL.value:
            expr = f"reshape(matmul({matrix(ins[0])},{matrix(ins[1])}),[{size}])"
        elif op in {Op.SCALE_ROWS.value, Op.SCALE_COLS.value}:
            diag, mat = ins if op == Op.SCALE_ROWS.value else (ins[1], ins[0])
            ds = offsets[diag][0]
            ms = offsets[mat][0]
            ts = offsets[i][0]
            b = graph.nodes[diag].shape[0]
            index = "ir" if op == Op.SCALE_ROWS.value else "ic"
            code += [
                f"do ic=0,{c - 1}",
                f"do ir=0,{r - 1}",
                f"w({ts}+ir+ic*{r})=w({ds}+{index}*(1+{b}))*w({ms}+ir+ic*{r})",
                "end do",
                "end do",
            ]
            expr = None
        elif op == Op.TRACE.value:
            b = graph.nodes[ins[0]].shape[0]
            a = offsets[ins[0]][0]
            ts = offsets[i][0]
            code += [f"w({ts})={_number(0)}", f"do it=0,{b - 1}", f"w({ts})=w({ts})+w({a}+it*(1+{b}))", "end do"]
            expr = None
        elif op == Op.NEG.value:
            expr = f"-{flat(ins[0])}"
        elif op == Op.CONJ.value:
            expr = f"conjg({flat(ins[0])})"
        elif op == Op.TRANSPOSE.value:
            expr = f"reshape(transpose({matrix(ins[0])}),[{size}])"
        elif op in {Op.NORMALIZE.value, Op.LOG_NORM.value}:
            code += [
                f"nm=maxval(abs({flat(ins[0])}))",
                "if (.not.ieee_is_finite(nm) .or. nm<=0.d0) then",
                "status=1",
                "return",
                "end if",
            ]
            expr = f"{flat(ins[0])}/nm" if op == Op.NORMALIZE.value else "cmplx(log(nm),0.d0,kind=c_double)"
        else:
            raise ValueError(f"unsupported operation {op}")
        if expr is not None:
            code.append(f"{dest}={expr}")
        # Preserve the two ordered multipliers just as the interpreter does.
        if node.scale != 1:
            code.append(f"{dest}={dest}*{_number(node.scale)}")
        if node.phase != 1:
            code.append(f"{dest}={dest}*{_number(node.phase)}")
    for name, i in graph.outputs.items():
        o = outputs[name]
        start = o["offset"] + 1
        end = o["offset"] + o["size"]
        code.append(f"y({start}:{end})={flat(i)}")
    code += [
        "if (.not.all(ieee_is_finite(real(y))) .or. .not.all(ieee_is_finite(aimag(y)))) status=2",
        "end subroutine",
        "end module",
    ]
    return Emitted("\n".join(code) + "\n", variables, outputs, n_input, n_work, n_output, len(live))


class CompiledGraph:
    def __init__(self, graph: Graph, directory: Path, compiler=None):
        self.graph = graph
        self.ir = emit_fortran(graph)
        compiler = compiler or shutil.which("gfortran")
        if not compiler:
            raise FileNotFoundError("gfortran unavailable")
        self.compiler = compiler
        directory = Path(directory)
        directory.mkdir(parents=True, exist_ok=True)
        digest = hashlib.sha256((self.ir.source + str(compiler) + " -O2 -fno-fast-math").encode()).hexdigest()[:20]
        self.source_path = directory / f"graph_{digest}.f90"
        self.library_path = directory / f"graph_{digest}.so"
        self.source_path.write_text(self.ir.source)
        start = time.perf_counter()
        # Each experiment directory is a reproducible build namespace, not a gate.
        proc = subprocess.run(
            [
                compiler,
                "-shared",
                "-fPIC",
                "-O2",
                "-fno-fast-math",
                "-ffree-line-length-none",
                "-J",
                str(directory),
                str(self.source_path),
                "-o",
                str(self.library_path),
            ],
            text=True,
            capture_output=True,
            timeout=180,
        )
        self.compile_s = time.perf_counter() - start
        if proc.returncode:
            raise RuntimeError(proc.stderr)
        self.lib = ctypes.CDLL(str(self.library_path.resolve()))
        self.fn = self.lib.scheme_d_eval
        self.fn.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.POINTER(ctypes.c_int)]
        self.fn.restype = None
        self.work = np.empty(max(1, self.ir.n_work), dtype=np.complex128)
        self.output = np.empty(self.ir.n_output, dtype=np.complex128)
        self.status = ctypes.c_int()

    def pack(self, env):
        for name in self.graph.unit_inputs:
            if name not in env or not np.array_equal(np.asarray(env[name]), np.ones((1, 1))):
                raise ValueError(f"unit-propagator guard failed: {name}")
        x = np.empty(max(1, self.ir.n_input), dtype=np.complex128)
        for name, item in self.ir.variables.items():
            shape, layout, *_ = item["spec"]
            value = _as_matrix(env[name], shape)
            if not np.isfinite(value).all():
                raise ValueError(f"nonfinite input {name}")
            if layout == Layout.DIAGONAL.value and not np.array_equal(value, np.diag(np.diag(value))):
                raise ValueError(f"diagonal input contract failed: {name}")
            x[item["offset"] : item["offset"] + item["size"]] = value.ravel(order="F")
        return x

    def run_packed(self, x):
        if x.dtype != np.complex128 or not x.flags.c_contiguous or x.size != max(1, self.ir.n_input):
            raise ValueError("invalid packed input buffer")
        self.fn(x.ctypes.data, self.work.ctypes.data, self.output.ctypes.data, ctypes.byref(self.status))
        if self.status.value:
            raise FloatingPointError(f"generated evaluator status={self.status.value}")
        return self.output

    def evaluate(self, env):
        result = self.run_packed(self.pack(env))
        return {
            name: result[o["offset"] : o["offset"] + o["size"]].reshape(o["shape"], order="F").copy()
            for name, o in self.ir.outputs.items()
        }
