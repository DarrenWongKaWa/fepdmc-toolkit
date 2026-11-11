"""Public C0 adapter contract (no network, no compiler)."""

from __future__ import annotations

import hashlib
import json
import os
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
APPLY = ROOT / "integration" / "c0" / "apply_c0.py"
VERIFY = ROOT / "integration" / "c0" / "verify_c0.py"
UPSTREAM_PIN = ROOT / "provenance" / "UPSTREAM_FEP_DMC.json"
C0_RECORD = ROOT / "provenance" / "C0_PUBLIC_PATCH.json"

TARGET = "perturbo-fep-dmc/pert-src/diagMC_JJ_updates.f90"
OTHER_PINNED = (
    "perturbo-fep-dmc/pert-src/diagMC.f90",
    "perturbo-fep-dmc/pert-src/pert_param.f90",
    "perturbo-fep-dmc/pert-src/makefile",
)
ANCHOR = "      call sample_q_omp_int(diagram%seed,vn1%i_q,vn1%Pq,vn1); Pq = vn1%Pq"
INSERTED = "      call cal_wq_int( vn1%i_q,     vn1%wq   )"
INTERNAL_SAMPLER = "      call sample_q_omp_int(diagram%seed, vn1%i_q, vn1%Pq, vn1); Pq = vn1%Pq"

PRISTINE_FORTRAN = """subroutine add_ph(diagram, stat)
      call sample_q_omp_int(diagram%seed, vn1%i_q, vn1%Pq, vn1); Pq = vn1%Pq
      call cal_wq_int( vn1%i_q,     vn1%wq   )
end subroutine

subroutine add_external_ph(diagram, stat)
      call sample_q_omp_int(diagram%seed,vn1%i_q,vn1%Pq,vn1); Pq = vn1%Pq

      vn2%wq = vn1%wq
end subroutine
"""


def _sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _write_layout(root: Path, payloads: dict[str, bytes]) -> None:
    for rel, data in payloads.items():
        path = root / rel
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(data)


def _init_git(root: Path) -> str:
    env = os.environ.copy()
    env["GIT_AUTHOR_NAME"] = "fepdmc-toolkit test"
    env["GIT_AUTHOR_EMAIL"] = "test@example.invalid"
    env["GIT_COMMITTER_NAME"] = env["GIT_AUTHOR_NAME"]
    env["GIT_COMMITTER_EMAIL"] = env["GIT_AUTHOR_EMAIL"]
    subprocess.check_call(["git", "init", "-q"], cwd=root, env=env)
    subprocess.check_call(["git", "config", "core.fileMode", "true"], cwd=root, env=env)
    subprocess.check_call(["git", "config", "core.autocrlf", "false"], cwd=root, env=env)
    subprocess.check_call(["git", "add", "-A"], cwd=root, env=env)
    subprocess.check_call(["git", "commit", "-q", "-m", "seed"], cwd=root, env=env)
    return subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=root, text=True).strip()


def _payloads(fortran: str, extra: dict[str, bytes] | None = None) -> dict[str, bytes]:
    out = {TARGET: fortran.encode("utf-8")}
    for rel in OTHER_PINNED:
        out[rel] = f"{rel}\n".encode()
    out["unrelated.txt"] = b"leave me alone\n"
    if extra:
        out.update(extra)
    return out


def _patched(fortran: str) -> str:
    assert fortran.count(ANCHOR + "\n") == 1
    return fortran.replace(ANCHOR + "\n", ANCHOR + "\n" + INSERTED + "\n", 1)


def _pin(payloads: dict[str, bytes], commit: str) -> dict:
    files = {rel: {"sha256": _sha(payloads[rel]), "role": "test"} for rel in (TARGET, *OTHER_PINNED)}
    return {
        "schema_version": 1,
        "role": "prospective_public_integration_baseline",
        "repository": "https://github.com/yaoluo/FEP-DMC",
        "reference_commit": commit,
        "reference_commit_date": "1970-01-01T00:00:00+00:00",
        "selected_reason": "unit-test fixture",
        "historical_donor_equivalence": "not_established",
        "files": files,
    }


def _c0_record(pre: bytes, post: bytes) -> dict:
    return {
        "schema_version": 1,
        "adapter": "c0_wq_refresh",
        "upstream_repository": "https://github.com/yaoluo/FEP-DMC",
        "upstream_commit": "test",
        "preimage_file": TARGET,
        "preimage_sha256": _sha(pre),
        "postimage_sha256": _sha(post),
        "changed_region": {
            "subroutine": "add_external_ph",
            "anchor": ANCHOR,
            "inserted_line": INSERTED,
        },
        "semantic_description": "test fixture",
    }


def _write_records(
    tmp: Path,
    payloads: dict[str, bytes],
    commit: str,
    post: bytes | None = None,
) -> tuple[Path, Path]:
    pre = payloads[TARGET]
    if post is None:
        post = _patched(pre.decode("utf-8")).encode("utf-8")
    pin = tmp / "pin.json"
    rec = tmp / "c0.json"
    pin.write_text(json.dumps(_pin(payloads, commit)))
    rec.write_text(json.dumps(_c0_record(pre, post)))
    return pin, rec


def _run_apply(tree: Path, pin: Path, rec: Path, *extra: str) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [sys.executable, str(APPLY), "--tree", str(tree), "--pin", str(pin), "--record", str(rec), *extra],
        capture_output=True,
        text=True,
    )


def _run_verify(tree: Path, pin: Path, rec: Path) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [sys.executable, str(VERIFY), "--tree", str(tree), "--pin", str(pin), "--record", str(rec)],
        capture_output=True,
        text=True,
    )


def test_c0_tools_and_record_exist():
    assert APPLY.is_file(), "integration/c0/apply_c0.py is required"
    assert VERIFY.is_file(), "integration/c0/verify_c0.py is required"
    assert C0_RECORD.is_file(), "provenance/C0_PUBLIC_PATCH.json is required"


def test_c0_record_matches_upstream_pin_and_denies_donor_equivalence():
    pin = json.loads(UPSTREAM_PIN.read_text())
    rec = json.loads(C0_RECORD.read_text())
    path = rec["preimage_file"]
    assert path == TARGET
    assert rec["upstream_repository"] == pin["repository"]
    assert rec["upstream_commit"] == pin["reference_commit"]
    assert rec["preimage_sha256"] == pin["files"][path]["sha256"]
    assert rec["preimage_sha256"] != rec["postimage_sha256"]
    assert len(rec["postimage_sha256"]) == 64
    assert rec["changed_region"]["anchor"] == ANCHOR
    assert rec["changed_region"]["inserted_line"] == INSERTED
    assert rec.get("historical_donor_equivalence", "not_established") != "established"


def test_dry_run_does_not_write(tmp_path: Path):
    payloads = _payloads(PRISTINE_FORTRAN)
    _write_layout(tmp_path, payloads)
    commit = _init_git(tmp_path)
    pin, rec = _write_records(tmp_path, payloads, commit)
    before = {rel: (tmp_path / rel).read_bytes() for rel in payloads}
    mtimes = {rel: (tmp_path / rel).stat().st_mtime_ns for rel in payloads}
    proc = _run_apply(tmp_path, pin, rec, "--dry-run")
    assert proc.returncode == 0, proc.stdout + proc.stderr
    assert "DRY_RUN" in proc.stdout
    assert INSERTED in proc.stdout
    assert (tmp_path / TARGET).read_bytes() == before[TARGET]
    for rel, data in before.items():
        assert (tmp_path / rel).read_bytes() == data
        assert (tmp_path / rel).stat().st_mtime_ns == mtimes[rel]
    leftover = list(tmp_path.rglob("*.fepdmc-toolkit-c0.tmp"))
    assert leftover == []


def test_apply_preserves_executable_mode(tmp_path: Path):
    payloads = _payloads(PRISTINE_FORTRAN)
    _write_layout(tmp_path, payloads)
    target = tmp_path / TARGET
    target.chmod(target.stat().st_mode | 0o111)
    commit = _init_git(tmp_path)
    pin, rec = _write_records(tmp_path, payloads, commit)
    proc = _run_apply(tmp_path, pin, rec)
    assert proc.returncode == 0, proc.stdout + proc.stderr
    assert target.stat().st_mode & 0o111


def test_apply_inserts_exactly_one_refresh(tmp_path: Path):
    payloads = _payloads(PRISTINE_FORTRAN)
    _write_layout(tmp_path, payloads)
    commit = _init_git(tmp_path)
    pin, rec = _write_records(tmp_path, payloads, commit)
    proc = _run_apply(tmp_path, pin, rec)
    assert proc.returncode == 0, proc.stdout + proc.stderr
    text = (tmp_path / TARGET).read_text()
    assert text.count(INSERTED) == 2  # add_ph already has the native call
    assert text == _patched(PRISTINE_FORTRAN)
    assert INTERNAL_SAMPLER in text
    rec_obj = json.loads(rec.read_text())
    assert _sha((tmp_path / TARGET).read_bytes()) == rec_obj["postimage_sha256"]
    for rel in OTHER_PINNED:
        assert (tmp_path / rel).read_bytes() == payloads[rel]
    assert (tmp_path / "unrelated.txt").read_bytes() == b"leave me alone\n"


def test_second_apply_is_idempotent(tmp_path: Path):
    payloads = _payloads(PRISTINE_FORTRAN)
    _write_layout(tmp_path, payloads)
    commit = _init_git(tmp_path)
    pin, rec = _write_records(tmp_path, payloads, commit)
    first = _run_apply(tmp_path, pin, rec)
    assert first.returncode == 0, first.stdout + first.stderr
    after_first = (tmp_path / TARGET).read_bytes()
    second = _run_apply(tmp_path, pin, rec)
    assert second.returncode == 2, second.stdout + second.stderr
    assert "ALREADY_APPLIED" in second.stdout
    assert (tmp_path / TARGET).read_bytes() == after_first
    assert after_first.decode().count(INSERTED) == 2


def test_wrong_preimage_hash_fails(tmp_path: Path):
    payloads = _payloads(PRISTINE_FORTRAN)
    _write_layout(tmp_path, payloads)
    commit = _init_git(tmp_path)
    pin, rec = _write_records(tmp_path, payloads, commit)
    obj = json.loads(rec.read_text())
    obj["preimage_sha256"] = "0" * 64
    rec.write_text(json.dumps(obj))
    proc = _run_apply(tmp_path, pin, rec)
    assert proc.returncode == 1
    assert (tmp_path / TARGET).read_bytes() == payloads[TARGET]


def test_missing_anchor_fails(tmp_path: Path):
    fortran = PRISTINE_FORTRAN.replace(ANCHOR, "      call not_the_sampler()")
    payloads = _payloads(fortran)
    _write_layout(tmp_path, payloads)
    commit = _init_git(tmp_path)
    dummy_post = payloads[TARGET] + b"!not-the-postimage\n"
    pin, rec = _write_records(tmp_path, payloads, commit, post=dummy_post)
    proc = _run_apply(tmp_path, pin, rec)
    assert proc.returncode == 1
    assert (tmp_path / TARGET).read_bytes() == payloads[TARGET]


def test_duplicate_anchor_fails(tmp_path: Path):
    fortran = PRISTINE_FORTRAN.replace(
        ANCHOR + "\n",
        ANCHOR + "\n" + ANCHOR + "\n",
        1,
    )
    payloads = _payloads(fortran)
    _write_layout(tmp_path, payloads)
    commit = _init_git(tmp_path)
    dummy_post = payloads[TARGET] + b"!not-the-postimage\n"
    pin, rec = _write_records(tmp_path, payloads, commit, post=dummy_post)
    proc = _run_apply(tmp_path, pin, rec)
    assert proc.returncode == 1
    assert (tmp_path / TARGET).read_bytes() == payloads[TARGET]


def test_partially_modified_source_fails(tmp_path: Path):
    payloads = _payloads(PRISTINE_FORTRAN)
    _write_layout(tmp_path, payloads)
    commit = _init_git(tmp_path)
    pin, rec = _write_records(tmp_path, payloads, commit)
    target = tmp_path / TARGET
    target.write_bytes(target.read_bytes() + b"! unexpected\n")
    proc = _run_apply(tmp_path, pin, rec)
    assert proc.returncode == 1
    v = _run_verify(tmp_path, pin, rec)
    assert v.returncode == 1
    assert "C0_STATE=UNKNOWN" in v.stdout


def test_verify_pristine_and_applied_and_unknown(tmp_path: Path):
    payloads = _payloads(PRISTINE_FORTRAN)
    _write_layout(tmp_path, payloads)
    commit = _init_git(tmp_path)
    pin, rec = _write_records(tmp_path, payloads, commit)
    pristine = _run_verify(tmp_path, pin, rec)
    assert pristine.returncode == 0, pristine.stdout + pristine.stderr
    assert "C0_STATE=PRISTINE" in pristine.stdout
    applied = _run_apply(tmp_path, pin, rec)
    assert applied.returncode == 0
    ok = _run_verify(tmp_path, pin, rec)
    assert ok.returncode == 0, ok.stdout + ok.stderr
    assert "C0_STATE=APPLIED" in ok.stdout
    (tmp_path / TARGET).write_bytes((tmp_path / TARGET).read_bytes() + b"x")
    bad = _run_verify(tmp_path, pin, rec)
    assert bad.returncode == 1
    assert "C0_STATE=UNKNOWN" in bad.stdout


def test_wrong_git_commit_is_refused(tmp_path: Path):
    payloads = _payloads(PRISTINE_FORTRAN)
    _write_layout(tmp_path, payloads)
    _init_git(tmp_path)
    pin, rec = _write_records(tmp_path, payloads, "b" * 40)
    proc = _run_apply(tmp_path, pin, rec)
    assert proc.returncode == 1
    assert (tmp_path / TARGET).read_bytes() == payloads[TARGET]


def test_git_unavailable_snapshot_is_refused(tmp_path: Path):
    payloads = _payloads(PRISTINE_FORTRAN)
    snapshot = tmp_path / "snapshot"
    _write_layout(snapshot, payloads)
    pin, rec = _write_records(tmp_path, payloads, "c" * 40)
    proc = _run_apply(snapshot, pin, rec)
    assert proc.returncode == 1
    assert (snapshot / TARGET).read_bytes() == payloads[TARGET]


def test_public_adapter_has_no_private_paths():
    public = APPLY.read_text()
    assert "night1" not in public
    assert "repair_v2" not in public
    lib = (ROOT / "integration" / "c0" / "c0_lib.py").read_text()
    assert "night1" not in lib
    assert "/Users/" not in public
    assert "/Users/" not in lib


def test_pristine_upstream_verifier_rejects_applied_tree(tmp_path: Path):
    payloads = _payloads(PRISTINE_FORTRAN)
    _write_layout(tmp_path, payloads)
    commit = _init_git(tmp_path)
    pin, rec = _write_records(tmp_path, payloads, commit)
    assert _run_apply(tmp_path, pin, rec).returncode == 0
    upstream = ROOT / "tools" / "verify_fep_dmc_upstream.py"
    proc = subprocess.run(
        [sys.executable, str(upstream), str(tmp_path), "--pin", str(pin)],
        capture_output=True,
        text=True,
    )
    assert proc.returncode != 0
    assert "FAIL" in proc.stdout
