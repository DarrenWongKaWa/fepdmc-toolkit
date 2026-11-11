"""Report PUBLIC_* P1/C0 composition state for a FEP-DMC checkout."""

from __future__ import annotations

import argparse

import p1_lib as p1


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    p1.add_tree_args(parser)
    args = parser.parse_args(argv)
    tree = args.tree.resolve()
    pin_path = args.pin.resolve()
    rec_path = args.record.resolve()
    c0_path = args.c0_record.resolve()
    if not tree.is_dir() or not pin_path.is_file() or not rec_path.is_file():
        print("P1_STATE=UNKNOWN")
        print("FAIL  tree/pin/record missing")
        return p1.EXIT_FAIL
    pin = p1.load_json(pin_path)
    rec = p1.load_json(rec_path)
    c0 = p1.load_json(c0_path) if c0_path.is_file() else None
    state, detail = p1.classify(tree, pin, rec, c0)
    print(f"P1_STATE={state}")
    print(detail)
    if state == p1.STATE_UNKNOWN:
        return p1.EXIT_FAIL
    if state == p1.STATE_PRISTINE:
        up = p1.run_upstream_verifier(tree, pin_path)
        if up.returncode != 0:
            print("FAIL  classified PRISTINE but upstream verifier failed")
            return p1.EXIT_FAIL
    else:
        print("pristine upstream verifier is expected to FAIL on this tree")
    return p1.EXIT_OK


if __name__ == "__main__":
    raise SystemExit(main())
