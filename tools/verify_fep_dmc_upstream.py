"""Fail-closed check of a public FEP-DMC tree against the pinned upstream commit.

Does not fetch, checkout, or modify the target.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
from pathlib import Path

EXIT_OK = 0
EXIT_FAIL = 1
EXIT_GIT_UNAVAILABLE = 3

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_PIN = ROOT / "provenance" / "UPSTREAM_FEP_DMC.json"


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _git(tree: Path, *args: str) -> subprocess.CompletedProcess[str]:
    argv = ["git", "-C", str(tree), *args]
    try:
        return subprocess.run(argv, capture_output=True, text=True)
    except FileNotFoundError:
        return subprocess.CompletedProcess(
            args=argv,
            returncode=127,
            stdout="",
            stderr="git executable not found",
        )


def _own_git_identity(tree: Path) -> bool:
    """True only when *tree* is a git work-tree toplevel.

    A source snapshot sitting inside some other repository must not
    inherit that repository's HEAD as FEP-DMC identity.
    """
    probe = _git(tree, "rev-parse", "--is-inside-work-tree")
    if probe.returncode != 0 or probe.stdout.strip() != "true":
        return False
    top = _git(tree, "rev-parse", "--show-toplevel")
    if top.returncode != 0 or not top.stdout.strip():
        return False
    try:
        return Path(top.stdout.strip()).resolve() == tree.resolve()
    except OSError:
        return False


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("tree", type=Path, help="FEP-DMC checkout or source snapshot")
    parser.add_argument(
        "--pin",
        type=Path,
        default=DEFAULT_PIN,
        help="pin JSON (default: provenance/UPSTREAM_FEP_DMC.json)",
    )
    args = parser.parse_args(argv)
    tree = args.tree.resolve()
    pin_path = args.pin.resolve()
    failed = False
    git_unavailable = False
    content_ok = True

    if not tree.is_dir():
        print(f"FAIL  tree missing: {tree}")
        return EXIT_FAIL
    if not pin_path.is_file():
        print(f"FAIL  pin missing: {pin_path}")
        return EXIT_FAIL

    rec = json.loads(pin_path.read_text())
    want_commit = rec["reference_commit"]
    files = rec["files"]

    have_git = _own_git_identity(tree)
    if have_git:
        got = _git(tree, "rev-parse", "HEAD")
        got_commit = got.stdout.strip()
        if got.returncode != 0 or len(got_commit) != 40:
            print("FAIL  git commit unreadable")
            failed = True
        elif got_commit != want_commit:
            print(f"FAIL  commit: got {got_commit} want {want_commit}")
            failed = True
        else:
            print(f"PASS  commit {want_commit}")
    else:
        git_unavailable = True
        print("GIT_UNAVAILABLE  no git identity in this tree")

    for rel, meta in files.items():
        path = tree / rel
        want = meta["sha256"]
        if not path.is_file():
            print(f"FAIL  missing {rel}")
            failed = True
            content_ok = False
            continue
        got_hash = _sha256(path)
        if got_hash != want:
            print(f"FAIL  {rel} sha256 mismatch")
            failed = True
            content_ok = False
            continue
        if have_git:
            st = _git(tree, "status", "--porcelain", "--", rel)
            if st.returncode != 0 or st.stdout.strip():
                print(f"FAIL  {rel} dirty working tree")
                failed = True
                continue
        print(f"PASS  {rel}")

    if git_unavailable and content_ok and not failed:
        print("CONTENT_MATCH  file hashes match; Git identity unavailable")
        return EXIT_GIT_UNAVAILABLE
    if failed:
        return EXIT_FAIL
    print("PASS  public FEP-DMC identity")
    return EXIT_OK


if __name__ == "__main__":
    raise SystemExit(main())
