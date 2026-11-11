"""Analyses of native FEP-DMC run directories.

Each analysis reads run directories written by "fepdmc-toolkit run" and prints
JSON. Most of them need the trace files of the research profile.

  fepdmc-toolkit analyze --list
  fepdmc-toolkit analyze <analysis> [arguments]
  fepdmc-toolkit analyze <analysis> --help

Every module also runs as "python -m fepdmc_toolkit.fepdmc.analysis.<analysis>"
and exposes its computation as a function that returns the result as a dict.
"""

from __future__ import annotations

import sys
from importlib import import_module

ANALYSES = (
    "analyze_mode_rb",
    "analyze_modes",
    "analyze_q",
    "analyze_rb",
    "analyze_svd",
    "collect_sign",
    "compare_bchain",
    "compare_between",
    "decompose_native",
    "estimate_groups",
    "summarize_materials",
    "summarize_sign",
)


def summary(name: str) -> str:
    """First line of the documentation of an analysis."""
    doc = import_module(f"{__name__}.{name}").__doc__ or ""
    return doc.strip().splitlines()[0] if doc.strip() else ""


def main(argv: list[str] | None = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    if not argv or argv[0] in ("-h", "--help"):
        print(__doc__)
        return 0 if argv else 2
    if argv[0] in ("-l", "--list"):
        for name in ANALYSES:
            print(f"{name:20s} {summary(name)}")
        return 0
    name = argv[0].replace("-", "_")
    if name not in ANALYSES:
        print(f"fepdmc-toolkit analyze: unknown analysis {argv[0]!r}; see --list", file=sys.stderr)
        return 2
    return int(import_module(f"{__name__}.{name}").main(argv[1:]) or 0)
