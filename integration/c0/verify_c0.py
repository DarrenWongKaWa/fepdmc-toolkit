"""Report C0_STATE=PRISTINE|APPLIED|UNKNOWN for a FEP-DMC checkout.

Does not fetch, checkout, or modify the target. Does not weaken
tools/verify_fep_dmc_upstream.py: a C0-patched tree is not a pristine
upstream identity.
"""

from __future__ import annotations

import argparse
import sys

import c0_lib as c0


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    c0.add_tree_args(parser)
    args = parser.parse_args(argv)
    tree = args.tree.resolve()
    pin_path = args.pin.resolve()
    rec_path = args.record.resolve()

    if not tree.is_dir():
        print(f"FAIL  tree missing: {tree}")
        print("C0_STATE=UNKNOWN")
        return c0.EXIT_FAIL
    if not pin_path.is_file() or not rec_path.is_file():
        print("FAIL  pin or C0 record missing")
        print("C0_STATE=UNKNOWN")
        return c0.EXIT_FAIL

    pin = c0.load_json(pin_path)
    rec = c0.load_json(rec_path)
    p1_path = args.p1_record.resolve()
    p1 = c0.load_json(p1_path) if p1_path.is_file() else None
    state, detail = c0.classify(tree, pin, rec, p1)
    print(f"C0_STATE={state}")
    print(detail)
    rel = rec["preimage_file"]
    path = tree / rel
    if path.is_file():
        print(f"file_sha256 {c0.sha256_path(path)}")
    print(f"preimage_sha256 {rec['preimage_sha256']}")
    print(f"postimage_sha256 {rec['postimage_sha256']}")

    if state == c0.STATE_PRISTINE:
        upstream = c0.run_upstream_verifier(tree, pin_path)
        sys.stdout.write(upstream.stdout)
        if upstream.returncode != 0:
            print("C0_STATE=UNKNOWN")
            print("FAIL  classify PRISTINE but upstream verifier did not pass")
            return c0.EXIT_FAIL
        return c0.EXIT_OK
    if state in (c0.STATE_APPLIED, c0.STATE_P1, c0.STATE_C0_P1):
        print("pristine upstream verifier is expected to FAIL on this tree")
        return c0.EXIT_OK
    return c0.EXIT_FAIL


if __name__ == "__main__":
    raise SystemExit(main())
