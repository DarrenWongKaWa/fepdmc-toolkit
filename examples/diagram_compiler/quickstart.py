"""Minimal installed-package example with a deterministic scalar binding."""

from __future__ import annotations

from fepdmc_toolkit.diagram_compiler import (
    Binding,
    build_scalar_ir,
    compile_cheap,
    compile_exact,
)


def main() -> None:
    ir = build_scalar_ir(1)
    x = Binding(k=0.0, q=(1.5707963267948966,), tau=(0.05, 0.12), t=1.0, omega=0.8, g=0.5, L=4)
    exact = compile_exact(ir).evaluate(x)["F"]
    cheap = compile_cheap(ir).evaluate(x)["F_hat"]
    print(f"family={ir.family.family_id}")
    print(f"exact_F={exact.real:.16g}")
    print(f"cheap_F_hat={cheap.real:.16g}")


if __name__ == "__main__":
    main()
