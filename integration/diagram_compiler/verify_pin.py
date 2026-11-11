"""Verify a pinned public FEP-DMC tree is PUBLIC_PRISTINE. Writes nothing."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import dc_lib as dc


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("tree", type=Path)
    parser.add_argument("--pin", type=Path, default=dc.DEFAULT_PIN)
    args = parser.parse_args(argv)
    tree = args.tree.resolve()
    pin = dc.load_json(args.pin.resolve())
    state, detail = dc.classify(tree, pin)
    print(f"DC_FEP_DMC_STATE={state}")
    print(detail)
    hits = dc.compiler_kernel_needles_in_tree(tree)
    if hits:
        print("FAIL  compiler-kernel needles in FEP-DMC tree:")
        for hit in hits:
            print(f"  {hit}")
        return dc.EXIT_FAIL
    print("NO_COMPILER_KERNEL_IN_FEP_DMC")
    if state != dc.STATE_PUBLIC_PRISTINE:
        return dc.EXIT_FAIL
    return dc.EXIT_OK


if __name__ == "__main__":
    sys.exit(main())
