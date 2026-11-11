# Validation records

| path | contents |
|---|---|
| `SHA256SUMS` | checksums of the validated artifacts (patched and generated Fortran, recorded results, upstream pins, these records); checked by `tests/test_integrity.py` |
| `diagram_compiler/` | reference artifacts of the diagram compiler: cheap policy, NativeKernelIR and generated Fortran of the two-band n = 1 kernel, `golden.json`, and the stage-1 and stage-6 differential results ([docs/diagram-compiler.md](../docs/diagram-compiler.md)) |
| `p1/` | native build and fixture validation of the P1 adapter ([docs/p1.md](../docs/p1.md)) |
| `fepdmc/` | records of the native FEP-DMC runs behind [docs/fepdmc.md](../docs/fepdmc.md) |

The records are kept as produced. Paths and names inside them refer to the
environment of the original runs.
