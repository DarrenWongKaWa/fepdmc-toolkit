# Validation stages of the diagram compiler

Each script repeats the checks of one stage
([docs/diagram-compiler.md](../docs/diagram-compiler.md)), writes its records
and `summary.json` to `--output` (default `build/validation/<stage>/`),
prints PASS or FAIL and exits with status 1 if a check fails.

| script | checks | needs |
|---|---|---|
| `dc_stage1_exact.py` | 224 bindings against the R1 reference and the index oracle | — |
| `dc_stage2_cheap.py` | cheap evaluator, antisymmetric score, recorded policy | — |
| `dc_stage3_da_kernel.py` | the three kernel paths; detailed balance on a three-state chain | — |
| `dc_stage4_native_ir.py` | NativeKernelIR equals the recorded IR; no exact operation on the stage-1 reject path | — |
| `dc_stage5_fortran_backend.py` | compiled kernel and primitives against the interpreter; emitted source equals the validated one | gfortran |
| `dc_stage6_fepdmc_toolchain.py` | the same kernel in the Ubuntu 22.04 image, next to a pristine clone of the pinned FEP-DMC commit | Docker, the clone ([integration/diagram_compiler/README.md](../integration/diagram_compiler/README.md)) |

```bash
for stage in scripts/dc_stage[1-5]_*.py; do python "$stage" || break; done
```
