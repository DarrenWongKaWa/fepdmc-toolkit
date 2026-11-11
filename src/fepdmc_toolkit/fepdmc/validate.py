"""End-to-end exactness check of the native fixes and the toolkit moves.

Runs independent kernels on the same target at low order, where every kernel
mixes well, and requires their pooled Q to agree:
  native  fixed native (FUTUREB_WQFIX, EXTFIX, EXTRMFIX), no toolkit move
  pure    only toolkit moves (general-span add/remove, ordered external
          add/remove, change-q on internal and external lines); native keeps
          phonon-mode changes only
  mixed   fixed native plus all toolkit moves
A kernel with a detailed-balance error, or one that samples a different
diagram space, shows up as a disagreement. The fully fixed build reproduced
the total energy E = -0.8763..-0.8767 at LiF-hole T = 500 K, maxOrder 7; without
the two removal fixes it gives -0.87311.

Usage:
  fepdmc-toolkit validate [--build TREE] [--chains 8] [--maxorder 7] [--tag val]
                           [--kernels native,pure,mixed] [--analyze-only]
                           [--ws DIR] [--data DIR] [--image IMAGE]
Each kernel is launched through fepdmc_toolkit.fepdmc.runner (workspace, dataset and
image as configured there). Exit status 0 if all pairs agree within --nsigma,
1 if a pair disagrees or a kernel has fewer than two finished chains, 2 if a run
could not be started.
"""

from __future__ import annotations

import argparse
import itertools
import math
import re
import subprocess
import sys
from pathlib import Path

import numpy as np

from . import runner

FIXES = ["--wqfix", "--extfix", "--extrmfix"]
KERNELS = {
    "native": FIXES + ["--chq", "1e-9"],
    "pure": ["--wqfix", "--pa", "0 0 1 0 0 0 0", "--any", "0.3", "--extar", "0.3", "--chq", "0.3", "--chq-ext"],
    "mixed": FIXES + ["--any", "0.1", "--extar", "0.1", "--chq", "0.05", "--chq-ext"],
}


def pooled(run_dir: Path, mat: str):
    """Pooled ratio estimate of E over chains and its chain-level delta-method SE."""
    E, S = [], []
    for c in sorted((run_dir / mat).glob("chain*")):
        s = (c / "stdout.log").read_text() if (c / "stdout.log").exists() else ""
        e = re.search(r"Etrue =\s+(\S+)", s)
        g = re.search(r"<g/Z> =\s+(\S+)", s)
        if e and g:
            E.append(float(e.group(1)))
            S.append(float(g.group(1)))
    if len(E) < 2:
        return float("nan"), float("nan"), len(E)
    E, S = np.array(E), np.array(S)
    n = E * S
    R = n.sum() / S.sum()
    h = (n - R * S) / S.mean()
    return float(R), float(h.std(ddof=1) / math.sqrt(len(E))), len(E)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        prog="fepdmc-toolkit validate", description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--build", default=runner.DEFAULT_BUILD, help="build tree name or container path of perturbo.x")
    runner.add_site_arguments(ap)
    ap.add_argument("--material", default="lif_hole", choices=list(runner.MATERIALS))
    ap.add_argument("--temp", type=float, default=500.0, help="temperature in K (default 500)")
    ap.add_argument("--maxorder", type=int, default=7, help="order cap (default 7)")
    ap.add_argument("--chains", type=int, default=8, help="chains per kernel (default 8)")
    ap.add_argument("--parallel", type=int, default=3, help="chains run at the same time per kernel (default 3)")
    ap.add_argument("--tag", default="validate", help="runs are written to <ws>/<tag>_<kernel>/")
    ap.add_argument("--kernels", default="native,pure,mixed", help="comma-separated subset of native,pure,mixed")
    ap.add_argument("--nsigma", type=float, default=3.0, help="largest accepted pairwise difference (default 3)")
    ap.add_argument("--analyze-only", action="store_true", help="compare existing runs; start no container")
    a = ap.parse_args(argv)
    names = [k.strip() for k in a.kernels.split(",") if k.strip()]
    unknown = [k for k in names if k not in KERNELS]
    if unknown or len(names) < 2:
        print(f"--kernels needs at least two of {', '.join(KERNELS)}; got {a.kernels!r}", file=sys.stderr)
        return 2
    if not a.analyze_only and not a.data:
        print("dataset directory required: --data or FEPDMC_TOOLKIT_DATA", file=sys.stderr)
        return 2
    procs = []
    if not a.analyze_only:
        for k in names:
            cmd = [
                sys.executable,
                "-m",
                "fepdmc_toolkit.fepdmc.runner",
                "ez",
                a.material,
                str(a.chains),
                "--build",
                a.build,
                "--ws",
                a.ws,
                "--image",
                a.image,
                "--temp",
                str(a.temp),
                "--maxorder",
                str(a.maxorder),
                "--tag",
                f"{a.tag}_{k}",
                "--parallel",
                str(a.parallel),
            ] + KERNELS[k]
            if a.data:
                cmd += ["--data", a.data]
            print("+", " ".join(cmd[1:]))
            procs.append(subprocess.Popen(cmd))
        codes = [p.wait() for p in procs]
        if any(codes):
            print("run failures:", codes, file=sys.stderr)
            return 2
    res = {k: pooled(Path(a.ws) / f"{a.tag}_{k}", a.material) for k in names}
    for k, (q, se, n) in res.items():
        print(f"{k:7s} n={n:2d}  E = {q:.5f} +- {se:.5f}")
    short = [k for k, (_, _, n) in res.items() if n < 2]
    if short:
        print(f"fewer than two finished chains for: {', '.join(short)}", file=sys.stderr)
        print("FAIL")
        return 1
    ok = True
    for x, y in itertools.combinations(names, 2):
        (qx, sx, _), (qy, sy, _) = res[x], res[y]
        z = (qx - qy) / math.hypot(sx, sy)
        good = abs(z) < a.nsigma
        ok &= good
        print(f"{x} vs {y}: {z:+.2f} sigma  {'ok' if good else 'DISAGREE'}")
    print("PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
