"""The compiler still emits exactly the artifacts that were validated natively.

The two-band n=1 delayed-acceptance kernel was compiled with gfortran in the
FEP-DMC toolchain image and checked against the Python oracle (stage 6). These
tests pin its NativeKernelIR, the emitted Fortran and the primitive runtime.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

from fepdmc_toolkit.diagram_compiler import (
    build_twoband_ir,
    compile_delayed_acceptance,
    emit_fortran,
    lower_native_kernel,
)
from fepdmc_toolkit.diagram_compiler.fortran_codegen import native_ir_digest
from fepdmc_toolkit.diagram_compiler.fortran_compile import RUNTIME_FILE

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "validation" / "diagram_compiler"
GOLDEN = json.loads((DATA / "golden.json").read_text())


def _kernel():
    return lower_native_kernel(compile_delayed_acceptance(build_twoband_ir(1)))


def test_emitted_fortran_is_the_validated_source():
    source = emit_fortran(_kernel())
    assert source == (DATA / "generated_twoband_n1.f90").read_text(encoding="utf-8")
    assert hashlib.sha256(source.encode("utf-8")).hexdigest() == GOLDEN["fortran_source_sha256"]


def test_native_ir_is_the_recorded_ir():
    nir = _kernel()
    assert native_ir_digest(nir) == GOLDEN["native_ir_digest"]
    assert json.loads(nir.to_json()) == json.loads((DATA / "NATIVE_KERNEL_IR.json").read_text())


def test_runtime_is_the_validated_runtime():
    assert hashlib.sha256(Path(RUNTIME_FILE).read_bytes()).hexdigest() == GOLDEN["runtime_sha256"]
