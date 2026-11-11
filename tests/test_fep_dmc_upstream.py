"""Public FEP-DMC upstream identity contract (no network)."""

from __future__ import annotations

import hashlib
import json
import os
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PIN = ROOT / "provenance" / "UPSTREAM_FEP_DMC.json"
TOOL = ROOT / "tools" / "verify_fep_dmc_upstream.py"

REQUIRED_FIELDS = (
    "schema_version",
    "role",
    "repository",
    "reference_commit",
    "reference_commit_date",
    "selected_reason",
    "historical_donor_equivalence",
    "files",
)

MUST_PIN_PATHS = (
    "perturbo-fep-dmc/pert-src/diagMC_JJ_updates.f90",
    "perturbo-fep-dmc/pert-src/diagMC.f90",
    "perturbo-fep-dmc/pert-src/pert_param.f90",
    "perturbo-fep-dmc/pert-src/makefile",
)


def _sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def test_pin_file_exists_with_required_fields():
    assert PIN.is_file(), "provenance/UPSTREAM_FEP_DMC.json is required"
    rec = json.loads(PIN.read_text())
    for key in REQUIRED_FIELDS:
        assert key in rec, key
    assert rec["repository"] == "https://github.com/yaoluo/FEP-DMC"
    assert rec["historical_donor_equivalence"] == "not_established"
    assert rec["role"] != "historical_donor"
    commit = rec["reference_commit"]
    assert isinstance(commit, str) and len(commit) == 40
    assert all(c in "0123456789abcdef" for c in commit)
    files = rec["files"]
    for path in MUST_PIN_PATHS:
        assert path in files, path
        digest = files[path]["sha256"]
        assert len(digest) == 64
        assert all(c in "0123456789abcdef" for c in digest)


def test_pin_does_not_claim_timed_donor_equivalence():
    rec = json.loads(PIN.read_text())
    assert rec["historical_donor_equivalence"] == "not_established"
    blob = json.dumps(rec).lower()
    assert "byte-identical to the private" not in blob


def test_verifier_script_exists():
    assert TOOL.is_file()


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


def _pin_for(payloads: dict[str, bytes], commit: str | None) -> dict:
    files = {rel: {"sha256": _sha(data), "role": "test"} for rel, data in payloads.items()}
    return {
        "schema_version": 1,
        "role": "prospective_public_integration_baseline",
        "repository": "https://github.com/yaoluo/FEP-DMC",
        "reference_commit": commit or ("a" * 40),
        "reference_commit_date": "1970-01-01T00:00:00+00:00",
        "selected_reason": "unit-test fixture",
        "historical_donor_equivalence": "not_established",
        "files": files,
    }


def _run(tree: Path, pin: Path) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [sys.executable, str(TOOL), str(tree), "--pin", str(pin)],
        capture_output=True,
        text=True,
    )


def test_verifier_passes_clean_matching_git_tree(tmp_path: Path):
    payloads = {p: f"{p}\n".encode() for p in MUST_PIN_PATHS}
    _write_layout(tmp_path, payloads)
    commit = _init_git(tmp_path)
    pin = tmp_path / "pin.json"
    pin.write_text(json.dumps(_pin_for(payloads, commit)))
    before = {rel: (tmp_path / rel).read_bytes() for rel in payloads}
    proc = _run(tmp_path, pin)
    assert proc.returncode == 0, proc.stdout + proc.stderr
    assert "PASS" in proc.stdout
    assert "FAIL" not in proc.stdout
    after = {rel: (tmp_path / rel).read_bytes() for rel in payloads}
    assert after == before


def test_verifier_fails_on_wrong_file_hash(tmp_path: Path):
    payloads = {p: f"{p}\n".encode() for p in MUST_PIN_PATHS}
    _write_layout(tmp_path, payloads)
    commit = _init_git(tmp_path)
    pin_obj = _pin_for(payloads, commit)
    pin_obj["files"][MUST_PIN_PATHS[0]]["sha256"] = "0" * 64
    pin = tmp_path / "pin.json"
    pin.write_text(json.dumps(pin_obj))
    proc = _run(tmp_path, pin)
    assert proc.returncode != 0
    assert "FAIL" in proc.stdout


def test_verifier_fails_on_missing_file(tmp_path: Path):
    payloads = {p: f"{p}\n".encode() for p in MUST_PIN_PATHS}
    missing = MUST_PIN_PATHS[0]
    keep = {k: v for k, v in payloads.items() if k != missing}
    _write_layout(tmp_path, keep)
    commit = _init_git(tmp_path)
    pin = tmp_path / "pin.json"
    pin.write_text(json.dumps(_pin_for(payloads, commit)))
    proc = _run(tmp_path, pin)
    assert proc.returncode != 0
    assert missing in proc.stdout


def test_verifier_fails_on_wrong_git_commit_even_if_hashes_match(tmp_path: Path):
    payloads = {p: f"{p}\n".encode() for p in MUST_PIN_PATHS}
    _write_layout(tmp_path, payloads)
    _init_git(tmp_path)
    pin = tmp_path / "pin.json"
    pin.write_text(json.dumps(_pin_for(payloads, "b" * 40)))
    proc = _run(tmp_path, pin)
    assert proc.returncode != 0
    assert "commit" in proc.stdout.lower()


def test_verifier_fails_on_dirty_pinned_file(tmp_path: Path):
    payloads = {p: f"{p}\n".encode() for p in MUST_PIN_PATHS}
    _write_layout(tmp_path, payloads)
    commit = _init_git(tmp_path)
    pin = tmp_path / "pin.json"
    pin.write_text(json.dumps(_pin_for(payloads, commit)))
    target = tmp_path / MUST_PIN_PATHS[0]
    target.write_bytes(target.read_bytes() + b"dirty\n")
    proc = _run(tmp_path, pin)
    assert proc.returncode == 1
    assert "FAIL" in proc.stdout


def test_verifier_fails_on_mode_dirty_file_even_if_hash_matches(tmp_path: Path):
    payloads = {p: f"{p}\n".encode() for p in MUST_PIN_PATHS}
    _write_layout(tmp_path, payloads)
    commit = _init_git(tmp_path)
    pin = tmp_path / "pin.json"
    pin.write_text(json.dumps(_pin_for(payloads, commit)))
    target = tmp_path / MUST_PIN_PATHS[0]
    target.chmod(target.stat().st_mode | 0o111)
    proc = _run(tmp_path, pin)
    assert proc.returncode == 1, proc.stdout + proc.stderr
    assert "dirty" in proc.stdout.lower()


def test_verifier_reports_git_unavailable_when_hashes_match(tmp_path: Path):
    payloads = {p: f"{p}\n".encode() for p in MUST_PIN_PATHS}
    _write_layout(tmp_path, payloads)
    pin = tmp_path / "pin.json"
    pin.write_text(json.dumps(_pin_for(payloads, "c" * 40)))
    proc = _run(tmp_path, pin)
    assert proc.returncode == 3
    assert "GIT_UNAVAILABLE" in proc.stdout
    assert "CONTENT_MATCH" in proc.stdout
    assert "FAIL" not in proc.stdout


def test_verifier_fails_closed_on_hash_mismatch_without_git(tmp_path: Path):
    payloads = {p: f"{p}\n".encode() for p in MUST_PIN_PATHS}
    _write_layout(tmp_path, payloads)
    pin_obj = _pin_for(payloads, "c" * 40)
    pin_obj["files"][MUST_PIN_PATHS[0]]["sha256"] = "0" * 64
    pin = tmp_path / "pin.json"
    pin.write_text(json.dumps(pin_obj))
    proc = _run(tmp_path, pin)
    assert proc.returncode == 1
    assert "FAIL" in proc.stdout
    assert "CONTENT_MATCH" not in proc.stdout


def test_nested_foreign_git_is_not_treated_as_fep_dmc_identity(tmp_path: Path):
    payloads = {p: f"{p}\n".encode() for p in MUST_PIN_PATHS}
    parent = tmp_path / "other_repo"
    parent.mkdir()
    (parent / "README").write_text("not fep-dmc\n")
    _init_git(parent)
    snapshot = parent / "snapshot"
    _write_layout(snapshot, payloads)
    pin = tmp_path / "pin.json"
    pin.write_text(json.dumps(_pin_for(payloads, "d" * 40)))
    proc = _run(snapshot, pin)
    assert proc.returncode == 3, proc.stdout + proc.stderr
    assert "GIT_UNAVAILABLE" in proc.stdout
    assert "CONTENT_MATCH" in proc.stdout
    assert "FAIL  commit" not in proc.stdout


def test_docs_quote_full_reference_commit():
    rec = json.loads(PIN.read_text())
    commit = rec["reference_commit"]
    assert commit in (ROOT / "docs" / "upstream.md").read_text()
    assert commit in (ROOT / "README.md").read_text()
    assert rec["historical_donor_equivalence"] == "not_established"
    assert rec["role"] == "prospective_public_integration_baseline"
