"""Collect the native average sign <g/Z> from existing FEP-DMC run directories.

Read-only over a local run tree (argument 1). For every directory holding
sign.dat-*, record the average sign sum_i sgn(Re D_i)/N written by
measure_gt (diagMC_gt.f90: gtrue = weight/abs(real(weight))), the number of
measurements Njjw, and the run settings needed to compare run groups.
Output paths are stored relative to the given root.
"""

from __future__ import annotations

import csv
import re
import sys
from pathlib import Path

from ._common import make_parser


def first_numbers(path: Path, line_no: int):
    lines = [l for l in path.read_text(errors="replace").splitlines() if l.strip()]
    return lines[line_no].split() if len(lines) > line_no else []


def namelist(path: Path, key: str):
    m = re.search(rf"^\s*{key}\s*=\s*([^!\n]+)", path.read_text(errors="replace"), re.M | re.I)
    return m.group(1).strip().strip("'\"") if m else ""


def collect(root: str, out: str) -> int:
    """Write the CSV ``out``; return the number of runs found under ``root``."""
    root_p = Path(root).resolve()
    rows = []
    for sign_file in sorted(root_p.rglob("sign.dat-*")):
        d = sign_file.parent
        vals = sign_file.read_text().split()
        row = {
            "run": str(d.relative_to(root_p)),
            "sign_re": float(vals[0]),
            "sign_im": float(vals[1]),
            "n_meas": int(vals[2]),
        }
        if (d / "diagMC.in").exists():
            nums = first_numbers(d / "diagMC.in", 1)
            row["Nmcmc_1e4"] = nums[0] if nums else ""
            row["maxOrder"] = nums[-1] if nums else ""
        if (d / "temper.in").exists():
            row["temper_first_line"] = " ".join(first_numbers(d / "temper.in", 0))
        if (d / "pert.in").exists():
            for k in ("prefix", "calc_mode", "DMC_Method", "band_min", "band_max"):
                row[k] = namelist(d / "pert.in", k)
        rows.append(row)
    fields = sorted({k for r in rows for k in r}, key=lambda k: (k != "run", k))
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)
    return len(rows)


def main(argv: list[str] | None = None) -> int:
    ap = make_parser(__name__, __doc__)
    ap.add_argument("root", help="directory tree with FEP-DMC run directories")
    ap.add_argument("out", help="CSV file to write")
    a = ap.parse_args(argv)
    if not Path(a.root).is_dir():
        print(f"error: {a.root}: not a directory", file=sys.stderr)
        return 2
    print(f"{collect(a.root, a.out)} runs -> {a.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
