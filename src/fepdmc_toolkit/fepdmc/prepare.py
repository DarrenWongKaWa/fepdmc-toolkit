"""Prepare a pristine FEP-DMC tree at the upstream pin in one step.

The ``base`` patch installs make.sys (gfortran/OpenMP/HDF5 flags for the QE 6.5
build) and the build hooks; the other patches follow in the order of
``fepdmc_toolkit.fepdmc.patches.ORDER``, each asserting its anchors exactly once.

Profiles:
  fixes     base + the three opt-in native fixes (wqfix, extfix, extrmfix)
  research  fixes + the toolkit moves and diagnostics (bchain), the default

The patches are first applied to a temporary copy of ``pert-src``. The tree is
written only if every patch applies, so a refused patch leaves it unchanged.
``--dry-run`` stops after the trial and lists the files that would be written.

Usage:
  fepdmc-toolkit prepare <perturbo-fep-dmc dir> [--profile fixes|research] [--dry-run]
"""

from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path

from . import FEP_DMC_PIN
from .patches import PROFILES, PatchError, apply

# Build products are not needed for the trial copy.
_SKIP = shutil.ignore_patterns("*.o", "*.mod", "*.a", "*.x")


@dataclass(frozen=True)
class Prepared:
    """Result of :func:`prepare_tree`."""

    profile: str
    patches: tuple[str, ...]
    files: tuple[str, ...]  # created or modified, relative to the perturbo-fep-dmc directory
    written: bool


def check_pin(root: Path) -> str:
    """Return a note on the tree's git revision against the upstream pin."""
    try:
        rev = subprocess.run(
            ["git", "-C", str(root), "rev-parse", "HEAD"], capture_output=True, text=True, check=True
        ).stdout.strip()
    except (subprocess.CalledProcessError, FileNotFoundError):
        return f"not a git checkout; expected pin {FEP_DMC_PIN}"
    return "pin ok" if rev == FEP_DMC_PIN else f"WARNING: tree at {rev}, expected pin {FEP_DMC_PIN}"


def _changed(stage: Path, root: Path) -> list[str]:
    """Files of ``stage`` that are new or differ from ``root``, as sorted relative paths."""
    out = []
    for path in sorted(stage.rglob("*")):
        if path.is_symlink() or not path.is_file():
            continue
        rel = path.relative_to(stage)
        target = root / rel
        if not target.is_file() or target.read_bytes() != path.read_bytes():
            out.append(rel.as_posix())
    return out


def prepare_tree(root: Path, profile: str = "research", dry_run: bool = False) -> Prepared:
    """Apply the patches of ``profile`` to the perturbo-fep-dmc directory ``root``.

    Raises :class:`PatchError` if a patch does not apply; ``root`` is then unchanged.
    With ``dry_run`` the patches are tried on a copy and nothing is written.
    """
    root = Path(root)
    if not (root / "pert-src").is_dir():
        raise FileNotFoundError(f"{root}: no pert-src/ directory")
    if profile not in PROFILES:
        raise ValueError(f"unknown profile {profile!r}; expected one of {sorted(PROFILES)}")
    steps = tuple(PROFILES[profile])
    with tempfile.TemporaryDirectory(prefix="fepdmc-toolkit-prepare-") as tmp:
        stage = Path(tmp) / "tree"
        shutil.copytree(root / "pert-src", stage / "pert-src", symlinks=True, ignore=_SKIP)
        if (root / "make.sys").is_file():
            shutil.copy2(root / "make.sys", stage / "make.sys")
        for name in steps:
            try:
                apply(stage, name)
            except PatchError as exc:
                raise PatchError(f"patch {name!r} does not apply: {exc}") from None
            except FileNotFoundError as exc:
                missing = Path(exc.filename).name if exc.filename else str(exc)
                raise PatchError(f"patch {name!r} does not apply: {missing} is missing") from None
        files = _changed(stage, root)
        if not dry_run:
            for rel in files:
                target = root / rel
                target.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy(stage / rel, target)  # content and permission bits
    return Prepared(profile=profile, patches=steps, files=tuple(files), written=not dry_run)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        prog="fepdmc-toolkit prepare", description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("root", type=Path, help="perturbo-fep-dmc directory (contains pert-src/)")
    ap.add_argument("--profile", choices=sorted(PROFILES), default="research", help="patch set (default: research)")
    ap.add_argument(
        "--dry-run", action="store_true", help="try the patches on a copy and list the files; write nothing"
    )
    a = ap.parse_args(argv)
    if not (a.root / "pert-src").is_dir():
        print(f"{a.root}: no pert-src/ directory", file=sys.stderr)
        return 2
    print(check_pin(a.root))
    try:
        result = prepare_tree(a.root, a.profile, a.dry_run)
    except PatchError as exc:
        print(f"{exc}", file=sys.stderr)
        print("No file was written. The patches apply only to an unpatched tree at the pinned commit.", file=sys.stderr)
        return 1
    would = "would " if a.dry_run else ""
    for name in result.patches:
        print(f"{would}apply {name}" if a.dry_run else f"applied {name}")
    for rel in result.files:
        print(f"{would}write {rel}" if a.dry_run else f"wrote {rel}")
    if a.dry_run:
        print(f"dry run ({a.profile}): the patches apply; {len(result.files)} files would be written.")
    else:
        print(
            f"done ({a.profile}): {len(result.files)} files written. Build pert-src with make inside the QE 6.5 tree."
        )
    return 0


if __name__ == "__main__":
    sys.exit(main())
