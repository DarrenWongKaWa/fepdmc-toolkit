"""Apply the public C0 wq-refresh to a pinned FEP-DMC checkout.

Does not fetch, checkout, or restore the target tree. Writes at most
the one pinned Fortran file, and only after the candidate postimage
hash matches the C0 record.
"""

from __future__ import annotations

import argparse
import sys

import c0_lib as c0


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    c0.add_tree_args(parser)
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="validate and print the diff; write nothing",
    )
    args = parser.parse_args(argv)
    tree = args.tree.resolve()
    pin_path = args.pin.resolve()
    rec_path = args.record.resolve()

    if not tree.is_dir():
        print(f"FAIL  tree missing: {tree}")
        return c0.EXIT_FAIL
    if not pin_path.is_file():
        print(f"FAIL  pin missing: {pin_path}")
        return c0.EXIT_FAIL
    if not rec_path.is_file():
        print(f"FAIL  C0 record missing: {rec_path}")
        return c0.EXIT_FAIL

    pin = c0.load_json(pin_path)
    rec = c0.load_json(rec_path)
    p1_path = args.p1_record.resolve()
    p1 = c0.load_json(p1_path) if p1_path.is_file() else None
    rel = rec["preimage_file"]
    path = tree / rel
    state, detail = c0.classify(tree, pin, rec, p1)
    print(f"C0_STATE={state}")

    if state in (c0.STATE_APPLIED, c0.STATE_C0_P1):
        print("ALREADY_APPLIED  exact C0 postimage; no write")
        return c0.EXIT_ALREADY
    if state not in (c0.STATE_PRISTINE, c0.STATE_P1):
        print(f"FAIL  refuse {state}: {detail}")
        return c0.EXIT_FAIL

    if state == c0.STATE_PRISTINE:
        upstream = c0.run_upstream_verifier(tree, pin_path)
        sys.stdout.write(upstream.stdout)
        if upstream.returncode != 0:
            print("FAIL  upstream identity is not official pristine")
            return c0.EXIT_FAIL

    source_bytes = path.read_bytes()
    if state == c0.STATE_PRISTINE and c0.sha256_bytes(source_bytes) != rec["preimage_sha256"]:
        print("FAIL  preimage SHA256 mismatch")
        return c0.EXIT_FAIL

    try:
        source = source_bytes.decode("utf-8")
    except UnicodeDecodeError:
        print("FAIL  target is not UTF-8")
        return c0.EXIT_FAIL

    try:
        after = c0.transform_source(
            source,
            anchor=rec["changed_region"]["anchor"],
            inserted_line=rec["changed_region"]["inserted_line"],
        )
    except c0.TransformError as exc:
        print(f"FAIL  {exc}")
        return c0.EXIT_FAIL

    after_bytes = after.encode("utf-8")
    post = c0.sha256_bytes(after_bytes)
    if state == c0.STATE_P1:
        want = (p1 or {}).get("states", {}).get("PUBLIC_C0_P1_APPLIED", {}).get(rel)
        if want is None or post != want:
            print(f"FAIL  C0-on-P1 candidate {post} != {want}")
            return c0.EXIT_FAIL
    elif post != rec["postimage_sha256"]:
        print(f"FAIL  candidate postimage {post} != record {rec['postimage_sha256']}")
        return c0.EXIT_FAIL

    print(c0.unified_diff(rel, source, after), end="")
    print(f"postimage_sha256 {post}")

    if args.dry_run:
        print("DRY_RUN  zero writes")
        return c0.EXIT_OK

    c0.atomic_write(path, after_bytes)
    if c0.sha256_path(path) != post:
        print("FAIL  write verification")
        return c0.EXIT_FAIL
    print(f"APPLIED  {rel}")
    print("revert with: git restore", rel)
    return c0.EXIT_OK


if __name__ == "__main__":
    raise SystemExit(main())
