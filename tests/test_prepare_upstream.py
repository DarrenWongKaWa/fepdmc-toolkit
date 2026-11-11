"""``fepdmc-toolkit prepare`` on a real clone of the pinned FEP-DMC commit.

Skipped unless FEPDMC_TOOLKIT_FEP_DMC names such a clone. The clone itself is not
modified: every test works on a copy of its perturbo-fep-dmc directory.
"""

from __future__ import annotations

import hashlib
import json
import os
import shutil
import subprocess
from pathlib import Path

import pytest

from fepdmc_toolkit.fepdmc import FEP_DMC_PIN, PatchError
from fepdmc_toolkit.fepdmc.prepare import check_pin, main, prepare_tree

ROOT = Path(__file__).resolve().parents[1]
RECORD = json.loads((ROOT / "provenance" / "PREPARED_TREE.json").read_text())
FEP = Path(os.environ["FEPDMC_TOOLKIT_FEP_DMC"]) if os.environ.get("FEPDMC_TOOLKIT_FEP_DMC") else None


def _at_pin(tree: Path | None) -> bool:
    if tree is None or not (tree / "perturbo-fep-dmc" / "pert-src").is_dir():
        return False
    head = subprocess.run(["git", "-C", str(tree), "rev-parse", "HEAD"], capture_output=True, text=True, check=False)
    return head.returncode == 0 and head.stdout.strip() == FEP_DMC_PIN


pytestmark = pytest.mark.skipif(not _at_pin(FEP), reason="clone of the pinned FEP-DMC commit not present")


def _snapshot(root: Path) -> dict[str, str]:
    return {
        p.relative_to(root).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest()
        for p in sorted(root.rglob("*"))
        if p.is_file() and not p.is_symlink()
    }


@pytest.fixture
def tree(tmp_path):
    dst = tmp_path / "perturbo-fep-dmc"
    shutil.copytree(FEP / "perturbo-fep-dmc", dst, symlinks=True)
    return dst


def test_record_names_the_pinned_commit():
    assert RECORD["reference_commit"] == FEP_DMC_PIN
    assert check_pin(FEP / "perturbo-fep-dmc") == "pin ok"
    assert set(RECORD["profiles"]["fixes"]) < set(RECORD["profiles"]["research"])


@pytest.mark.parametrize("profile", ["fixes", "research"])
def test_prepare_reproduces_the_recorded_tree(profile, tree):
    want = RECORD["profiles"][profile]
    before = _snapshot(tree)

    dry = prepare_tree(tree, profile, dry_run=True)
    assert not dry.written and sorted(dry.files) == sorted(want)
    assert _snapshot(tree) == before

    result = prepare_tree(tree, profile)
    assert result.written and result.files == dry.files
    after = _snapshot(tree)
    assert {rel: after[rel] for rel in want} == want
    assert {rel for rel in after if after[rel] != before.get(rel)} == set(want)
    assert set(before) <= set(after)

    with pytest.raises(PatchError):  # a patched tree is refused, and stays as it is
        prepare_tree(tree, profile)
    assert _snapshot(tree) == after


def test_prepare_command_reports_the_files(tree, capsys):
    assert main([str(tree), "--profile", "fixes", "--dry-run"]) == 0
    out = capsys.readouterr().out
    assert out.count("would write ") == len(RECORD["profiles"]["fixes"])
    assert main([str(tree), "--profile", "fixes"]) == 0
    assert "13 files written" in capsys.readouterr().out
    assert main([str(tree), "--profile", "research"]) == 1
    assert "No file was written" in capsys.readouterr().err
