"""Shared C0 source-transform helpers. Not a general patch engine."""

from __future__ import annotations

import argparse
import difflib
import hashlib
import importlib.util
import json
import os
import re
import subprocess
import sys
from pathlib import Path

EXIT_OK = 0
EXIT_FAIL = 1
EXIT_ALREADY = 2

STATE_PRISTINE = "PRISTINE"
STATE_APPLIED = "APPLIED"
STATE_UNKNOWN = "UNKNOWN"

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
DEFAULT_PIN = ROOT / "provenance" / "UPSTREAM_FEP_DMC.json"
DEFAULT_RECORD = ROOT / "provenance" / "C0_PUBLIC_PATCH.json"
DEFAULT_P1 = ROOT / "provenance" / "P1_PUBLIC_PATCH.json"
UPSTREAM_TOOL = ROOT / "tools" / "verify_fep_dmc_upstream.py"
JJ_UPDATES = "perturbo-fep-dmc/pert-src/diagMC_JJ_updates.f90"
MAKEFILE = "perturbo-fep-dmc/pert-src/makefile"
JJ_DRIVER = "perturbo-fep-dmc/pert-src/diagMC_JJ.f90"
MODULE_REL = "perturbo-fep-dmc/pert-src/linear_da_mod.f90"

STATE_P1 = "P1_APPLIED"
STATE_C0_P1 = "C0_P1_APPLIED"

SUBROUTINE_START = re.compile(r"^[ \t]*subroutine\s+add_external_ph\s*\(", re.M | re.I)
SUBROUTINE_END = re.compile(r"^[ \t]*end subroutine\b", re.M | re.I)
ANCHOR_RE = re.compile(
    r"^[ \t]*call sample_q_omp_int\(diagram%seed,vn1%i_q,vn1%Pq,vn1\); Pq = vn1%Pq[ \t]*\r?\n",
    re.M,
)


class TransformError(ValueError):
    pass


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha256_path(path: Path) -> str:
    return sha256_bytes(path.read_bytes())


def load_json(path: Path) -> dict:
    return json.loads(path.read_text())


def load_upstream_mod():
    spec = importlib.util.spec_from_file_location("verify_fep_dmc_upstream", UPSTREAM_TOOL)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load {UPSTREAM_TOOL}")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def run_upstream_verifier(tree: Path, pin: Path) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [sys.executable, str(UPSTREAM_TOOL), str(tree), "--pin", str(pin)],
        capture_output=True,
        text=True,
    )


def newline_for(source: str) -> str:
    return "\r\n" if "\r\n" in source else "\n"


def external_body_span(source: str) -> tuple[int, int]:
    starts = list(SUBROUTINE_START.finditer(source))
    if len(starts) != 1:
        raise TransformError(f"expected exactly one add_external_ph, found {len(starts)}")
    rest = source[starts[0].end() :]
    ends = list(SUBROUTINE_END.finditer(rest))
    if not ends:
        raise TransformError("missing end subroutine for add_external_ph")
    return starts[0].start(), starts[0].end() + ends[0].end()


def transform_source(source: str, *, anchor: str, inserted_line: str) -> str:
    start, end = external_body_span(source)
    body = source[start:end]
    matches = list(ANCHOR_RE.finditer(body))
    if len(matches) == 0:
        raise TransformError("missing add_external_ph sample_q_omp_int anchor")
    if len(matches) > 1:
        raise TransformError("duplicate add_external_ph sample_q_omp_int anchor")
    matched = matches[0].group(0).strip()
    if matched != anchor.strip():
        raise TransformError("anchor text does not match the C0 record")
    nl = newline_for(source)
    insertion = inserted_line + nl
    at = matches[0].end()
    remainder = body[at:]
    if remainder.startswith(insertion):
        raise TransformError("refresh call already present")
    if re.search(r"call\s+cal_wq_int\s*\(", body, re.I):
        raise TransformError("add_external_ph already contains cal_wq_int")
    new_body = body[:at] + insertion + remainder
    return source[:start] + new_body + source[end:]


def unified_diff(rel: str, before: str, after: str) -> str:
    return "".join(
        difflib.unified_diff(
            before.splitlines(keepends=True),
            after.splitlines(keepends=True),
            fromfile=f"a/{rel}",
            tofile=f"b/{rel}",
            n=3,
        )
    )


def _file_dirty(mod, tree: Path, rel: str) -> bool:
    st = mod._git(tree, "status", "--porcelain", "--", rel)
    return st.returncode != 0 or bool(st.stdout.strip())


def classify(tree: Path, pin: dict, rec: dict, p1: dict | None = None) -> tuple[str, str]:
    """Return (state, detail). Official C0 identity requires own git toplevel."""
    mod = load_upstream_mod()
    target_rel = rec["preimage_file"]
    pre = rec["preimage_sha256"]
    post = rec["postimage_sha256"]
    want_commit = pin["reference_commit"]
    files = pin["files"]
    if target_rel not in files:
        return STATE_UNKNOWN, f"pin missing {target_rel}"
    if files[target_rel]["sha256"] != pre:
        return STATE_UNKNOWN, "C0 preimage does not match upstream pin"

    have_git = mod._own_git_identity(tree)
    if not have_git:
        return STATE_UNKNOWN, "git identity unavailable"

    got = mod._git(tree, "rev-parse", "HEAD")
    got_commit = got.stdout.strip()
    if got.returncode != 0 or len(got_commit) != 40 or got_commit != want_commit:
        return STATE_UNKNOWN, "git commit is not the pinned upstream"

    hashes: dict[str, str] = {}
    for rel in files:
        path = tree / rel
        if not path.is_file():
            return STATE_UNKNOWN, f"missing {rel}"
        hashes[rel] = sha256_path(path)
    for rel in (JJ_DRIVER, MODULE_REL):
        extra = tree / rel
        if extra.is_file():
            hashes[rel] = sha256_path(extra)

    if p1 and "states" in p1:
        st = p1["states"]

        def match_p1_state(name: str) -> bool:
            wanted = st.get(name) or {}
            if not wanted:
                return False
            return all(hashes.get(rel) == digest for rel, digest in wanted.items())

        def pin_untouched_ok() -> str | None:
            skip = {JJ_UPDATES, MAKEFILE}
            for rel in files:
                if rel in skip:
                    continue
                if hashes[rel] != files[rel]["sha256"]:
                    return f"{rel} hash mismatch"
                if _file_dirty(mod, tree, rel):
                    return f"{rel} dirty working tree"
            return None

        if match_p1_state("PUBLIC_P1_APPLIED"):
            want_mod = p1.get("module_sha256")
            if want_mod and hashes.get(MODULE_REL) != want_mod:
                return STATE_UNKNOWN, "P1 module identity mismatch"
            why = pin_untouched_ok()
            if why:
                return STATE_UNKNOWN, why
            return STATE_P1, "P1 applied, C0 not applied"
        if match_p1_state("PUBLIC_C0_P1_APPLIED"):
            want_mod = p1.get("module_sha256")
            if want_mod and hashes.get(MODULE_REL) != want_mod:
                return STATE_UNKNOWN, "P1 module identity mismatch"
            why = pin_untouched_ok()
            if why:
                return STATE_UNKNOWN, why
            return STATE_C0_P1, "C0 and P1 applied"

    others = [rel for rel in files if rel != target_rel]
    for rel in others:
        if hashes[rel] != files[rel]["sha256"]:
            return STATE_UNKNOWN, f"{rel} hash mismatch"
        if _file_dirty(mod, tree, rel):
            return STATE_UNKNOWN, f"{rel} dirty working tree"

    target_hash = hashes[target_rel]
    if target_hash == pre:
        if _file_dirty(mod, tree, target_rel):
            return STATE_UNKNOWN, f"{target_rel} dirty working tree"
        return STATE_PRISTINE, "pinned pristine source"
    if target_hash == post:
        return STATE_APPLIED, "exact C0 postimage"
    return STATE_UNKNOWN, "target matches neither preimage nor postimage"


def atomic_write(path: Path, data: bytes) -> None:
    tmp = path.with_name(path.name + ".fepdmc-toolkit-c0.tmp")
    mode = path.stat().st_mode
    try:
        tmp.write_bytes(data)
        os.chmod(tmp, mode)
        os.replace(tmp, path)
    finally:
        if tmp.exists():
            try:
                tmp.unlink()
            except OSError:
                pass


def add_tree_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--tree", type=Path, required=True, help="FEP-DMC checkout")
    parser.add_argument("--pin", type=Path, default=DEFAULT_PIN)
    parser.add_argument("--record", type=Path, default=DEFAULT_RECORD)
    parser.add_argument("--p1-record", type=Path, default=DEFAULT_P1)
