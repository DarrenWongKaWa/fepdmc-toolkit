"""Average native sign per run group, from the table written by collect_sign.

Runs are grouped by the directory that holds them (the part of the run path
before ``/runs/``, else the parent directory). For each group and for all runs
together: number of chains, mean sign and its standard error from the spread of
independent chains. ``variance_ceiling_perfect_grouping`` is 1/<s>^2, the largest
variance gain a grouped measure could reach by removing every sign cancellation.

Runs with fewer than --min-nmcmc steps (in units of 1e4) are skipped, so that
short test runs do not enter the average.
"""

from __future__ import annotations

import csv
import math
import statistics
import sys
from collections import defaultdict
from pathlib import Path

from ._common import AnalysisError, emit, make_parser


def group_of(run: str) -> str:
    """Run group of a run path as written by collect_sign."""
    return run.split("/runs/")[0] if "/runs/" in run else run.rsplit("/", 1)[0]


def analyze(table: str, min_nmcmc: int = 100, exclude: tuple[str, ...] = ()) -> dict:
    with open(table, newline="") as fh:
        rows = list(csv.DictReader(fh))
    groups = defaultdict(list)
    for r in rows:
        g = group_of(r["run"])
        if any(pattern in g for pattern in exclude) or int(r.get("Nmcmc_1e4") or 0) < min_nmcmc:
            continue
        groups[g].append(float(r["sign_re"]))
    if not groups:
        raise AnalysisError(f"{table}: no run with at least {min_nmcmc}e4 steps")
    out = {"groups": {}}
    allv = []
    for g, v in sorted(groups.items()):
        se = statistics.stdev(v) / math.sqrt(len(v)) if len(v) > 1 else float("nan")
        out["groups"][g] = {"chains": len(v), "mean_sign": statistics.fmean(v), "stderr": se}
        allv += v
    s = statistics.fmean(allv)
    out["pooled"] = {
        "chains": len(allv),
        "mean_sign": s,
        "stderr": statistics.stdev(allv) / math.sqrt(len(allv)) if len(allv) > 1 else float("nan"),
    }
    out["variance_ceiling_perfect_grouping"] = 1.0 / s**2
    return out


def main(argv: list[str] | None = None) -> int:
    ap = make_parser(__name__, __doc__)
    ap.add_argument("table", help="CSV file written by collect_sign")
    ap.add_argument(
        "--min-nmcmc", type=int, default=100, help="skip runs with fewer steps, in units of 1e4 (default 100)"
    )
    ap.add_argument(
        "--exclude", action="append", default=[], metavar="TEXT", help="skip run groups whose path contains TEXT"
    )
    a = ap.parse_args(argv)
    if not Path(a.table).is_file():
        print(f"error: {a.table}: not a file", file=sys.stderr)
        return 2
    return emit(lambda: analyze(a.table, a.min_nmcmc, tuple(a.exclude)))


if __name__ == "__main__":
    sys.exit(main())
