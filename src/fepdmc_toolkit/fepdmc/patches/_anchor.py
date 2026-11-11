"""Anchored text edits shared by the source patches.

Every edit names the exact text it replaces. A patch that does not find its
anchor exactly once raises :class:`PatchError` instead of guessing, so it either
applies cleanly to the pinned, unpatched tree or does not apply at all.
"""

from __future__ import annotations

from pathlib import Path


class PatchError(RuntimeError):
    """The tree is not the pinned, unpatched source that the patch expects."""


def read(path: Path) -> str:
    """Read a source file as UTF-8, independent of the locale."""
    return Path(path).read_text(encoding="utf-8")


def write(path: Path, text: str) -> None:
    """Write a source file as UTF-8 with LF line endings."""
    Path(path).write_text(text, encoding="utf-8", newline="\n")


def locate(text: str, marker: str, start: int = 0) -> int:
    """Offset of ``marker`` in ``text`` at or after ``start``."""
    pos = text.find(marker, start)
    if pos < 0:
        raise PatchError(f"marker {marker!r} not found")
    return pos


def once(text: str, old: str, new: str, label: str) -> str:
    """Replace the single occurrence of ``old`` in ``text``."""
    n = text.count(old)
    if n != 1:
        raise PatchError(f"anchor for {label!r} found {n} times, expected exactly once")
    return text.replace(old, new)
