"""Validated artifacts are byte-identical to the versions that were validated.

validation/SHA256SUMS lists the Fortran sources that are patched into FEP-DMC or
compiled, the recorded benchmark results, the upstream pins and the validation
records. Editing any of them invalidates the native validation; update the
checksum only together with a new validation run.
"""

from __future__ import annotations

import hashlib
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SUMS = ROOT / "validation" / "SHA256SUMS"


def _entries():
    for line in SUMS.read_text().splitlines():
        if line.strip() and not line.startswith("#"):
            digest, rel = line.split(maxsplit=1)
            yield digest, rel.strip()


def test_checksum_file_covers_the_patched_and_generated_fortran():
    listed = {rel for _, rel in _entries()}
    for rel in (
        "src/fepdmc_toolkit/fepdmc/data/bchain.f90",
        "src/fepdmc_toolkit/fepdmc/data/rb_window.f90",
        "integration/p1/linear_da_mod.f90",
        "src/fepdmc_toolkit/diagram_compiler/fortran/native_kernel_runtime_v1.f90",
        "validation/diagram_compiler/generated_twoband_n1.f90",
        "benchmarks/p1_vs_bbest/frozen_results.json",
        "provenance/UPSTREAM_FEP_DMC.json",
    ):
        assert rel in listed, rel


def test_validated_artifacts_are_unchanged():
    changed = [rel for digest, rel in _entries() if hashlib.sha256((ROOT / rel).read_bytes()).hexdigest() != digest]
    assert changed == []
