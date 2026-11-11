"""Public P1 source transforms. Not a generic patch engine."""

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

STATE_PRISTINE = "PUBLIC_PRISTINE"
STATE_C0 = "PUBLIC_C0_APPLIED"
STATE_P1 = "PUBLIC_P1_APPLIED"
STATE_C0_P1 = "PUBLIC_C0_P1_APPLIED"
STATE_UNKNOWN = "UNKNOWN"

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
DEFAULT_PIN = ROOT / "provenance" / "UPSTREAM_FEP_DMC.json"
DEFAULT_C0 = ROOT / "provenance" / "C0_PUBLIC_PATCH.json"
DEFAULT_P1 = ROOT / "provenance" / "P1_PUBLIC_PATCH.json"
UPSTREAM_TOOL = ROOT / "tools" / "verify_fep_dmc_upstream.py"
CANONICAL_MODULE = HERE / "linear_da_mod.f90"

JJ_UPDATES = "perturbo-fep-dmc/pert-src/diagMC_JJ_updates.f90"
JJ_DRIVER = "perturbo-fep-dmc/pert-src/diagMC_JJ.f90"
MAKEFILE = "perturbo-fep-dmc/pert-src/makefile"
MODULE_REL = "perturbo-fep-dmc/pert-src/linear_da_mod.f90"

SWAP_START = re.compile(r"^[ \t]*subroutine\s+update_swap\s*\(", re.M | re.I)
SUB_END = re.compile(r"^[ \t]*end subroutine\b", re.M | re.I)
TIME_STOP_RE = re.compile(
    r"^[ \t]*if\(tauL>tauR\) stop ' time disordered @ update_swap'[ \t]*\r?\n",
    re.M,
)
MAKE_FRAG = "diagMC.f90 \\\ndiagMC_debug.f90 \\"
MAKE_FRAG_NEW = "diagMC.f90 \\\nlinear_da_mod.f90 \\\ndiagMC_debug.f90 \\"

USE_LINE = (
    "      use linear_da_mod, only: linear_da_is_on, linear_da_eval, &\n"
    "           linear_da_aux_uniform, linear_da_stage2, linear_da_note_expensive, &\n"
    "           linear_da_note_accept, linear_da_note_stage2_reject, linear_da_note_s1_reject\n"
)
DECL_LINE = "      logical :: da_el, da_acc\n      real(dp) :: ell_hat, u1\n"
STAGE1_BLOCK = """      da_el = .false.
      ell_hat = 0.d0
      if (linear_da_is_on()) then
         call linear_da_eval(diagram, iv1, iv2, da_el, ell_hat)
         if (.not. da_el) then
            stop 'linear_da: ineligible swap with LINEAR_DA=on'
         endif
         call linear_da_aux_uniform(u1)
         if (log(max(u1, 1.0e-300_dp)) >= min(0.d0, ell_hat)) then
            call linear_da_note_s1_reject()
            return
         endif
      endif
"""
NOTE_EXPENSIVE = "      if (linear_da_is_on()) call linear_da_note_expensive()\n"
STAGE2_BLOCK = """      if (linear_da_is_on()) then
         if (.not. (real(mat_old) == real(mat_old)) .or. abs(real(mat_old)) == 0.d0) then
            stop 'linear_da: Re(mat_old) is 0 or nonfinite'
         endif
         call linear_da_stage2(ran, P_accept, ell_hat, da_acc)
      else
         da_acc = (ran < P_accept)
      endif

      if(da_acc) then
"""
ACCEPT_NOTE = "         if (linear_da_is_on()) call linear_da_note_accept()\n"
S2_ELSE = "      else\n         if (linear_da_is_on()) call linear_da_note_stage2_reject()\n"


class TransformError(ValueError):
    pass


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha256_path(path: Path) -> str:
    return sha256_bytes(path.read_bytes())


def load_json(path: Path) -> dict:
    return json.loads(path.read_text())


def newline_for(source: str) -> str:
    return "\r\n" if "\r\n" in source else "\n"


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


def swap_span(source: str) -> tuple[int, int]:
    starts = list(SWAP_START.finditer(source))
    if len(starts) != 1:
        raise TransformError(f"expected exactly one update_swap, found {len(starts)}")
    rest = source[starts[0].end() :]
    ends = list(SUB_END.finditer(rest))
    if not ends:
        raise TransformError("missing end subroutine for update_swap")
    return starts[0].start(), starts[0].end() + ends[0].end()


def _nl(source: str, block: str) -> str:
    if newline_for(source) == "\r\n":
        return block.replace("\n", "\r\n")
    return block


def transform_updates(source: str) -> str:
    start, end = swap_span(source)
    body = source[start:end]
    if "linear_da_eval" in body:
        raise TransformError("update_swap already contains P1 hooks")
    stops = list(TIME_STOP_RE.finditer(body))
    if len(stops) == 0:
        raise TransformError("missing update_swap time-order stop")
    if len(stops) > 1:
        raise TransformError("duplicate update_swap time-order stop")
    use_pat = re.compile(r"^[ \t]*use DiagMC[ \t]*\r?\n", re.M)
    uses = list(use_pat.finditer(body))
    if len(uses) != 1:
        raise TransformError("expected exactly one use DiagMC in update_swap")
    body = body[: uses[0].end()] + _nl(source, USE_LINE) + body[uses[0].end() :]
    implicit = re.search(r"^[ \t]*implicit none[ \t]*\r?\n", body, re.M)
    if implicit is None:
        raise TransformError("missing implicit none in update_swap")
    body = body[: implicit.end()] + _nl(source, DECL_LINE) + body[implicit.end() :]
    stops = list(TIME_STOP_RE.finditer(body))
    body = body[: stops[0].end()] + _nl(source, STAGE1_BLOCK) + body[stops[0].end() :]
    gkq = list(re.finditer(r"^[ \t]*call cal_gkq_vtex_int\( vRn, vRn%gkq \)[ \t]*\r?\n", body, re.M))
    if len(gkq) != 1:
        raise TransformError("expected unique first cal_gkq_vtex_int(vRn) in update_swap")
    body = body[: gkq[0].start()] + _nl(source, NOTE_EXPENSIVE) + body[gkq[0].start() :]
    ran = list(
        re.finditer(
            r"^[ \t]*call random_number_omp\(diagram%seed,ran\)[ \t]*\r?\n"
            r"[ \t]*\r?\n"
            r"^[ \t]*if\(ran<P_accept\) then[ \t]*\r?\n",
            body,
            re.M,
        )
    )
    if len(ran) != 1:
        raise TransformError("expected unique native ran < P_accept block")
    prefix = re.match(
        r"([ \t]*call random_number_omp\(diagram%seed,ran\)[ \t]*\r?\n[ \t]*\r?\n)",
        ran[0].group(0),
    )
    if prefix is None:
        raise TransformError("could not split native accept block")
    body = body[: ran[0].start()] + prefix.group(1) + _nl(source, STAGE2_BLOCK) + body[ran[0].end() :]
    acc = list(re.finditer(r"^[ \t]*stat%accept\(7\) = stat%accept\(7\) \+ 1[ \t]*\r?\n", body, re.M))
    if len(acc) != 1:
        raise TransformError("expected unique stat%accept(7) increment in update_swap")
    body = body[: acc[0].end()] + _nl(source, ACCEPT_NOTE) + body[acc[0].end() :]
    term = list(
        re.finditer(
            r"^[ \t]*end if[ \t]*\r?\n(?:[ \t]*\r?\n)?[ \t]*end subroutine",
            body,
            re.M,
        )
    )
    if len(term) != 1:
        raise TransformError("expected unique terminal end if of update_swap accept")
    body = body[: term[0].start()] + _nl(source, S2_ELSE) + term[0].group(0)
    if "linear_da_stage2" not in body:
        raise TransformError("missing linear_da_stage2 call")
    if body.find("linear_da_eval") > body.find("cal_gkq_vtex_int"):
        raise TransformError("Stage 1 is after cal_gkq")
    return source[:start] + body + source[end:]


def transform_makefile(source: str) -> str:
    if "linear_da_mod.f90" in source:
        raise TransformError("makefile already lists linear_da_mod.f90")
    count = source.count(MAKE_FRAG)
    if count != 1:
        raise TransformError(f"makefile PERTMOD fragment count={count}")
    return source.replace(MAKE_FRAG, MAKE_FRAG_NEW, 1)


def transform_driver(source: str) -> str:
    if "linear_da_ensure" in source:
        raise TransformError("diagMC_JJ.f90 already contains P1 hooks")
    use = list(re.finditer(r"^[ \t]*use DiagMC[ \t]*\r?\n", source, re.M))
    if len(use) != 1:
        raise TransformError("expected exactly one use DiagMC in diagMC_JJ.f90")
    ins_use = "   use linear_da_mod, only: linear_da_ensure, linear_da_report\n"
    source = source[: use[0].end()] + _nl(source, ins_use) + source[use[0].end() :]
    setup = list(re.finditer(r"^[ \t]*call setup_dqmc\(\)[ \t]*\r?\n", source, re.M))
    if len(setup) != 1:
        raise TransformError("expected unique call setup_dqmc() in diagMC_JJ.f90")
    source = source[: setup[0].end()] + _nl(source, "   call linear_da_ensure()\n") + source[setup[0].end() :]
    acc = list(
        re.finditer(
            r"^[ \t]*write\(stdout,'\(A30,7E15\.5\)'\)'acceptance = ',.*\r?\n",
            source,
            re.M,
        )
    )
    if len(acc) != 1:
        raise TransformError("expected unique acceptance write in diagMC_JJ.f90")
    source = source[: acc[0].end()] + _nl(source, "   call linear_da_report()\n") + source[acc[0].end() :]
    if source.find("call setup_dqmc()") > source.find("call linear_da_ensure()"):
        raise TransformError("ensure is not after setup_dqmc")
    if source.find("'acceptance = '") > source.find("call linear_da_report()"):
        raise TransformError("report is not after acceptance write")
    return source


def _file_dirty(mod, tree: Path, rel: str) -> bool:
    st = mod._git(tree, "status", "--porcelain", "--", rel)
    return st.returncode != 0 or bool(st.stdout.strip())


def _hash_if(tree: Path, rel: str) -> str | None:
    path = tree / rel
    if not path.is_file():
        return None
    return sha256_path(path)


def classify(tree: Path, pin: dict, p1: dict, c0: dict | None) -> tuple[str, str]:
    mod = load_upstream_mod()
    want = pin["reference_commit"]
    if not mod._own_git_identity(tree):
        return STATE_UNKNOWN, "git identity unavailable"
    got = mod._git(tree, "rev-parse", "HEAD")
    commit = got.stdout.strip()
    if got.returncode != 0 or len(commit) != 40 or commit != want:
        return STATE_UNKNOWN, "git commit is not the pinned upstream"
    pin_files = pin["files"]
    hashes = {}
    for rel in pin_files:
        path = tree / rel
        if not path.is_file():
            return STATE_UNKNOWN, f"missing {rel}"
        hashes[rel] = sha256_path(path)
    for rel in (JJ_DRIVER, MODULE_REL):
        h = _hash_if(tree, rel)
        if h is not None:
            hashes[rel] = h
    states = p1["states"]

    def match(name: str) -> bool:
        st = states[name]
        for rel, digest in st.items():
            if hashes.get(rel) != digest:
                return False
        return True

    others = [rel for rel in pin_files if rel not in (JJ_UPDATES, MAKEFILE)]
    for rel in others:
        if hashes[rel] != pin_files[rel]["sha256"]:
            return STATE_UNKNOWN, f"{rel} hash mismatch"
        if _file_dirty(mod, tree, rel):
            return STATE_UNKNOWN, f"{rel} dirty working tree"

    if match(STATE_PRISTINE):
        for rel in (JJ_UPDATES, MAKEFILE):
            if _file_dirty(mod, tree, rel):
                return STATE_UNKNOWN, f"{rel} dirty working tree"
        if MODULE_REL in hashes:
            return STATE_UNKNOWN, "unexpected linear_da_mod.f90 on pristine tree"
        return STATE_PRISTINE, "vanilla pin"
    if match(STATE_C0):
        if MODULE_REL in hashes:
            return STATE_UNKNOWN, "unexpected module on C0-only tree"
        return STATE_C0, "C0 applied"
    if match(STATE_P1):
        if hashes.get(MODULE_REL) != p1["module_sha256"]:
            return STATE_UNKNOWN, "linear_da_mod.f90 is not the canonical module"
        return STATE_P1, "P1 applied"
    if match(STATE_C0_P1):
        if hashes.get(MODULE_REL) != p1["module_sha256"]:
            return STATE_UNKNOWN, "linear_da_mod.f90 is not the canonical module"
        return STATE_C0_P1, "C0 and P1 applied"
    return STATE_UNKNOWN, "hashes match no recorded public state"


def atomic_write(path: Path, data: bytes) -> None:
    tmp = path.with_name(path.name + ".fepdmc-toolkit-p1.tmp")
    mode = path.stat().st_mode if path.exists() else 0o644
    try:
        path.parent.mkdir(parents=True, exist_ok=True)
        tmp.write_bytes(data)
        os.chmod(tmp, mode)
        os.replace(tmp, path)
    finally:
        if tmp.exists():
            try:
                tmp.unlink()
            except OSError:
                pass


def commit_writes(items: list[tuple[Path, bytes]]) -> None:
    staged: list[tuple[Path, Path]] = []
    try:
        for path, data in items:
            tmp = path.with_name(path.name + ".fepdmc-toolkit-p1.tmp")
            path.parent.mkdir(parents=True, exist_ok=True)
            tmp.write_bytes(data)
            if path.exists():
                os.chmod(tmp, path.stat().st_mode)
            staged.append((tmp, path))
        for tmp, path in staged:
            os.replace(tmp, path)
    finally:
        for tmp, _path in staged:
            if tmp.exists():
                try:
                    tmp.unlink()
                except OSError:
                    pass


def add_tree_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--tree", type=Path, required=True)
    parser.add_argument("--pin", type=Path, default=DEFAULT_PIN)
    parser.add_argument("--record", type=Path, default=DEFAULT_P1)
    parser.add_argument("--c0-record", type=Path, default=DEFAULT_C0)
