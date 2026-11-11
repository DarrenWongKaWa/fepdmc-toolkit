"""Fail-closed FEP-DMC pin checks for the validation of the diagram compiler.

Does not modify the FEP-DMC tree. In stage 6 the validation executable is built
outside the FEP-DMC tree; the pinned tree is only checked to be pristine and
free of the compiler kernel.
"""

from __future__ import annotations

import hashlib
import json
import subprocess
from pathlib import Path

EXIT_OK = 0
EXIT_FAIL = 1

STATE_PUBLIC_PRISTINE = "PUBLIC_PRISTINE"
STATE_UNKNOWN = "UNKNOWN"

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
DEFAULT_PIN = ROOT / "provenance" / "UPSTREAM_FEP_DMC.json"
PAYLOAD_REL = "perturbo-fep-dmc/pert-src/dc_validate"


def sha256_path(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def git_head(tree: Path) -> str | None:
    probe = subprocess.run(
        ["git", "-C", str(tree), "rev-parse", "HEAD"],
        capture_output=True,
        text=True,
        check=False,
    )
    if probe.returncode != 0:
        return None
    return probe.stdout.strip()


def own_git_toplevel(tree: Path) -> bool:
    inside = subprocess.run(
        ["git", "-C", str(tree), "rev-parse", "--is-inside-work-tree"],
        capture_output=True,
        text=True,
        check=False,
    )
    if inside.returncode != 0 or inside.stdout.strip() != "true":
        return False
    top = subprocess.run(
        ["git", "-C", str(tree), "rev-parse", "--show-toplevel"],
        capture_output=True,
        text=True,
        check=False,
    )
    if top.returncode != 0:
        return False
    try:
        return Path(top.stdout.strip()).resolve() == tree.resolve()
    except OSError:
        return False


def production_hashes(tree: Path, pin: dict) -> dict[str, str]:
    out = {}
    for rel in pin["files"]:
        path = tree / rel
        out[rel] = sha256_path(path) if path.is_file() else ""
    return out


def classify(tree: Path, pin: dict) -> tuple[str, str]:
    if not tree.is_dir():
        return STATE_UNKNOWN, "tree missing"
    if not own_git_toplevel(tree):
        return STATE_UNKNOWN, "not an own git toplevel"
    head = git_head(tree)
    want = pin["reference_commit"]
    if head != want:
        return STATE_UNKNOWN, f"HEAD {head} != {want}"
    for rel, meta in pin["files"].items():
        path = tree / rel
        if not path.is_file():
            return STATE_UNKNOWN, f"missing {rel}"
        got = sha256_path(path)
        if got != meta["sha256"]:
            return STATE_UNKNOWN, f"{rel} hash mismatch"
    extra = tree / PAYLOAD_REL
    if extra.exists():
        return STATE_UNKNOWN, "unexpected validation payload directory in FEP-DMC tree"
    return STATE_PUBLIC_PRISTINE, "pin files and HEAD match"


def compiler_kernel_needles_in_tree(tree: Path) -> list[str]:
    hits = []
    pert = tree / "perturbo-fep-dmc" / "pert-src"
    if not pert.is_dir():
        return ["missing pert-src"]
    needles = (
        "evaluate_da_kernel",
        "native_kernel_runtime_v1",
        "generated_da_kernel",
        "dc_validate",
    )
    for path in pert.rglob("*"):
        if not path.is_file():
            continue
        if path.suffix.lower() not in {".f90", ".f", ".makefile", ""} and path.name != "makefile":
            continue
        try:
            text = path.read_text(encoding="utf-8", errors="ignore")
        except OSError:
            continue
        for needle in needles:
            if needle in text:
                hits.append(f"{path.relative_to(tree)}:{needle}")
    return hits
