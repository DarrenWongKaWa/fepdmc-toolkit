"""Apply public P1 delayed acceptance to a pinned FEP-DMC checkout.

Does not fetch, checkout, or restore the target. Writes the four
P1 files only after every candidate postimage matches the record.
"""

from __future__ import annotations

import argparse
import sys

import p1_lib as p1


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    p1.add_tree_args(parser)
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args(argv)
    tree = args.tree.resolve()
    pin_path = args.pin.resolve()
    rec_path = args.record.resolve()
    c0_path = args.c0_record.resolve()

    if not tree.is_dir():
        print(f"FAIL  tree missing: {tree}")
        return p1.EXIT_FAIL
    if not pin_path.is_file() or not rec_path.is_file():
        print("FAIL  pin or P1 record missing")
        return p1.EXIT_FAIL

    pin = p1.load_json(pin_path)
    rec = p1.load_json(rec_path)
    c0 = p1.load_json(c0_path) if c0_path.is_file() else None
    state, detail = p1.classify(tree, pin, rec, c0)
    print(f"P1_STATE={state}")

    if state in (p1.STATE_P1, p1.STATE_C0_P1):
        print("ALREADY_APPLIED  exact P1 postimage; no write")
        return p1.EXIT_ALREADY
    if state not in (p1.STATE_PRISTINE, p1.STATE_C0):
        print(f"FAIL  refuse {state}: {detail}")
        return p1.EXIT_FAIL

    if state == p1.STATE_PRISTINE:
        up = p1.run_upstream_verifier(tree, pin_path)
        sys.stdout.write(up.stdout)
        if up.returncode != 0:
            print("FAIL  upstream identity is not official pristine")
            return p1.EXIT_FAIL
        want_name = p1.STATE_P1
    else:
        want_name = p1.STATE_C0_P1

    try:
        upd = (tree / p1.JJ_UPDATES).read_text(encoding="utf-8")
        mk = (tree / p1.MAKEFILE).read_text(encoding="utf-8")
        drv = (tree / p1.JJ_DRIVER).read_text(encoding="utf-8")
        new_upd = p1.transform_updates(upd)
        new_mk = p1.transform_makefile(mk)
        new_drv = p1.transform_driver(drv)
    except (OSError, UnicodeDecodeError, p1.TransformError) as exc:
        print(f"FAIL  {exc}")
        return p1.EXIT_FAIL

    module = p1.CANONICAL_MODULE.read_bytes()
    if p1.sha256_bytes(module) != rec["module_sha256"]:
        print("FAIL  canonical linear_da_mod.f90 does not match the P1 record")
        return p1.EXIT_FAIL

    candidates = {
        p1.JJ_UPDATES: new_upd.encode("utf-8"),
        p1.MAKEFILE: new_mk.encode("utf-8"),
        p1.JJ_DRIVER: new_drv.encode("utf-8"),
        p1.MODULE_REL: module,
    }
    want = rec["states"][want_name]
    for rel, digest in want.items():
        got = p1.sha256_bytes(candidates[rel])
        if got != digest:
            print(f"FAIL  candidate {rel} {got} != {digest}")
            return p1.EXIT_FAIL
    if p1.sha256_bytes(module) != rec["module_sha256"]:
        print("FAIL  module hash")
        return p1.EXIT_FAIL

    print(p1.unified_diff(p1.JJ_UPDATES, upd, new_upd), end="")
    print(p1.unified_diff(p1.MAKEFILE, mk, new_mk), end="")
    print(p1.unified_diff(p1.JJ_DRIVER, drv, new_drv), end="")
    print(f"would add {p1.MODULE_REL}")
    for rel, data in candidates.items():
        print(f"postimage_sha256 {rel} {p1.sha256_bytes(data)}")

    if args.dry_run:
        print("DRY_RUN  zero writes")
        return p1.EXIT_OK

    writes = [(tree / rel, data) for rel, data in candidates.items()]
    p1.commit_writes(writes)
    for rel, data in candidates.items():
        if p1.sha256_path(tree / rel) != p1.sha256_bytes(data):
            print(f"FAIL  write verification {rel}")
            return p1.EXIT_FAIL
    print(f"APPLIED  {want_name}")
    print(
        "revert with: git restore",
        p1.JJ_UPDATES,
        p1.MAKEFILE,
        p1.JJ_DRIVER,
        ";",
        "rm",
        p1.MODULE_REL,
    )
    return p1.EXIT_OK


if __name__ == "__main__":
    raise SystemExit(main())
