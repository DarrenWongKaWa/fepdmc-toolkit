"""Native EZ average sign and formation energy per material.

For every <runs>/<material>/chain*/stdout.log read <g/Z> (the average sign
under DMC_Method 0), Etrue, the bare band energy printed as 'Ek-mu' (meV)
and wall time. Q = Etrue - E_bare. Errors are the spread of independent
chains (FUTUREB_SEED per chain). Also reports the sign ceiling 1/<s>^2 on
any grouped-measure variance gain.
"""

from __future__ import annotations

import json
import math
import re
import statistics
import sys
from pathlib import Path

from ._common import AnalysisError, emit, make_parser


def parse(stdout: str) -> dict:
    g = re.search(r"<g/Z> =\s+(\S+)", stdout)
    e = re.search(r"Etrue =\s+(\S+)", stdout)
    b = re.search(r"Ek-mu\s*=\s*(\S+)\(meV\)", stdout)
    beta = re.search(r"beta =\s+(\S+)", stdout)
    return {
        "sign": float(g.group(1)) if g else None,
        "Etrue": float(e.group(1)) if e else None,
        "E_bare": float(b.group(1)) / 1000.0 if b else None,
        "beta": float(beta.group(1)) if beta else None,
    }


def analyze(runs: str) -> dict:
    res = {}
    for mdir in sorted(Path(runs).iterdir()):
        chains = []
        for c in sorted(mdir.glob("chain*")):
            so = c / "stdout.log"
            if not so.exists():
                continue
            d = parse(so.read_text())
            t = re.search(r"real\s+(\d+)m([\d.]+)s", (c / "stderr.log").read_text())
            d["wall_s"] = int(t.group(1)) * 60 + float(t.group(2)) if t else None
            d["meta"] = json.loads((c / "run_meta.json").read_text())
            if d["sign"] is not None:
                chains.append(d)
        if not chains:
            continue
        s = [c["sign"] for c in chains]
        q = [c["Etrue"] - c["E_bare"] for c in chains]
        n = len(chains)
        res[mdir.name] = {
            "chains": n,
            "bands": chains[0]["meta"]["bands"],
            "hole": chains[0]["meta"]["hole"],
            "T_K": chains[0]["meta"]["T"],
            "beta_eV-1": chains[0]["beta"],
            "mean_sign": statistics.fmean(s),
            "se_sign": statistics.stdev(s) / math.sqrt(n) if n > 1 else None,
            "min_sign": min(s),
            "max_sign": max(s),
            "Q_eV": statistics.fmean(q),
            "se_Q_eV": statistics.stdev(q) / math.sqrt(n) if n > 1 else None,
            "variance_ceiling_perfect_grouping": 1.0 / statistics.fmean(s) ** 2,
            "mean_wall_s": statistics.fmean(c["wall_s"] for c in chains if c["wall_s"]),
        }
    if not res:
        raise AnalysisError(f"{runs}: no material directory with finished chains")
    return res


def main(argv: list[str] | None = None) -> int:
    ap = make_parser(__name__, __doc__)
    ap.add_argument("runs", help="run directory with one subdirectory per material")
    ap.add_argument("out", nargs="?", default=None, help="also write the JSON to this file")
    a = ap.parse_args(argv)

    def compute():
        res = analyze(a.runs)
        if a.out:
            Path(a.out).write_text(json.dumps(res, indent=2))
        return res

    return emit(compute)


if __name__ == "__main__":
    sys.exit(main())
