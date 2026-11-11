"""Packaging: version agreement, sdist contents, and no private paths in the runtime package."""

from __future__ import annotations

import re
import tomllib
from pathlib import Path

import fepdmc_toolkit

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src" / "fepdmc_toolkit"

# Directories the tests read from; the sdist must contain them so that pytest runs from an unpacked sdist.
SDIST_DIRS = (
    "benchmarks",
    "data",
    "docs",
    "examples",
    "integration",
    "provenance",
    "scripts",
    "tests",
    "tools",
    "validation",
)

PRIVATE_PATH_NEEDLES = ("research/", "CLOSURE", "night1", "C2_fixture", "/Users/")


def test_version_is_defined_once_and_cited_consistently():
    version = fepdmc_toolkit.__version__
    config = tomllib.loads((ROOT / "pyproject.toml").read_text())
    assert config["project"]["dynamic"] == ["version"]
    assert config["tool"]["setuptools"]["dynamic"]["version"] == {"attr": "fepdmc_toolkit.__version__"}
    citation = re.search(r'^version: "([^"]+)"', (ROOT / "CITATION.cff").read_text(), re.M).group(1)
    assert citation == version
    assert "## " + version in (ROOT / "CHANGELOG.md").read_text()


def test_manifest_grafts_the_test_inputs():
    body = (ROOT / "MANIFEST.in").read_text()
    missing = [d for d in SDIST_DIRS if f"graft {d}" not in body]
    assert missing == []


def test_runtime_python_has_no_private_paths():
    hits = []
    for path in SRC.rglob("*.py"):
        text = path.read_text()
        for needle in PRIVATE_PATH_NEEDLES:
            if needle in text:
                hits.append(f"{path.relative_to(ROOT)}: {needle}")
    assert hits == []


def test_console_script_points_at_the_cli():
    project = tomllib.loads((ROOT / "pyproject.toml").read_text())["project"]
    assert project["scripts"] == {"fepdmc-toolkit": "fepdmc_toolkit.fepdmc.cli:main"}
