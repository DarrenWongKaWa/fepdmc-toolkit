"""Offline estimate of wider mode groups (crossing lines carried) from lines_dump.dat.

For every dumped native sample: the internal lines (la, lb) and their exact
single-line mode ratio rho_l. A group rule chooses a line set S from topology
only; lines of S crossed by no other S line are folded (L), the others (C) carry
their mode index through a left-to-right contraction, recursing into folded
lines. Width = max over cuts of the C lines open at that cut within the same
innermost folded container. Rule "K": greedy by opening order, add a line if the
width stays <= K (K=0 is the non-crossing set).

Estimates, using the validated per-line independence (pair tests ~1e-3):
  Z_B/Z_A ~ E_A[prod_{l in S} rho_l],   <s>_B ~ <s>_A / E_A[...]
plus a cost proxy sum over positions of n^{width at that position}.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

from ._common import AnalysisError, emit, make_parser


def parse(path: Path):
    samples = []
    cur = None
    for line in path.read_text().splitlines():
        if line.startswith("#"):
            continue
        t = line.split()
        if t[0] == "M":
            cur = dict(sign=float(t[3]), lines=[])
            samples.append(cur)
        else:
            cur["lines"].append((int(t[0]), int(t[1]), float(t[2])))
    return samples


def crosses(x, y):
    (a, b), (c, d) = x, y
    return a < c < b < d or c < a < d < b


def classify(S, lines):
    L, C = [], []
    for i in S:
        if any(j != i and crosses(lines[i][:2], lines[j][:2]) for j in S):
            C.append(i)
        else:
            L.append(i)
    return L, C


def width_and_cost(S, lines, n):
    L, C = classify(S, lines)
    if not S:
        return 0, 0.0
    npos = max(lines[i][1] for i in S) + 1
    # innermost folded container of each cut x (between positions x and x+1)
    cont = [-1] * (npos + 1)
    best = [10**9] * (npos + 1)
    for i in L:
        a, b = lines[i][:2]
        for x in range(a, b):
            if b - a < best[x]:
                best[x] = b - a
                cont[x] = i

    # container of a C line = innermost L line containing it entirely
    def container_of(i):
        a, b = lines[i][:2]
        c, bs = -1, 10**9
        for j in L:
            p, q = lines[j][:2]
            if p < a and b < q and q - p < bs:
                c, bs = j, q - p
        return c

    ccont = {i: container_of(i) for i in C}
    w = 0
    cost = 0.0
    for x in range(npos):
        k = sum(1 for i in C if lines[i][0] <= x < lines[i][1] and ccont[i] == cont[x])
        w = max(w, k)
        cost += n**k
    return w, cost


def select(lines, K, n):
    S = []
    for i in range(len(lines)):
        trial = S + [i]
        wt, _ = width_and_cost(trial, lines, n)
        if wt <= K:
            S = trial
    return S


def analyze(runs: str, mat: str, n: int = 6) -> dict:
    Ks = [0, 1, 2, 3]
    out = {"material": mat, "n_modes": n, "rules": {}}
    per_chain = {}
    for f in sorted((Path(runs) / mat).glob("chain*/lines_dump.dat")):
        samples = parse(f)
        sA = np.mean([s["sign"] for s in samples])
        rec = {
            "n_samples": len(samples),
            "sign_A": float(sA),
            "mean_lines": float(np.mean([len(s["lines"]) for s in samples])),
        }
        for K in Ks:
            prods, fr, cost = [], [], []
            for s in samples:
                lines = s["lines"]
                S = select(lines, K, n)
                prods.append(float(np.prod([lines[i][2] for i in S])) if S else 1.0)
                fr.append(len(S) / max(1, len(lines)))
                cost.append(width_and_cost(S, lines, n)[1] / max(1, max((l[1] for l in lines), default=1)))
            e = float(np.mean(prods))
            rec[f"K{K}"] = dict(
                fraction_of_lines=float(np.mean(fr)),
                E_prod_rho=e,
                sign_B_est=float(sA / e),
                cost_per_position=float(np.mean(cost)),
            )
        e_all = float(np.mean([np.prod([l[2] for l in s["lines"]]) if s["lines"] else 1.0 for s in samples]))
        rec["all_lines"] = dict(E_prod_rho=e_all, sign_B_est=float(sA / e_all))
        per_chain[f.parent.name] = rec
    if not per_chain:
        raise AnalysisError(f"{Path(runs) / mat}: no chain with lines_dump.dat")
    out["per_chain"] = per_chain
    for key in [f"K{K}" for K in Ks] + ["all_lines"]:
        vals = [c[key] for c in per_chain.values()]
        sA = np.mean([c["sign_A"] for c in per_chain.values()])
        e = np.mean([v["E_prod_rho"] for v in vals])
        out["rules"][key] = dict(
            E_prod_rho=float(e),
            sign_B_est=float(sA / e),
            fraction_of_lines=float(np.mean([v.get("fraction_of_lines", 1.0) for v in vals])),
            cost_per_position=float(np.mean([v.get("cost_per_position", float("nan")) for v in vals])),
        )
    out["sign_A"] = float(np.mean([c["sign_A"] for c in per_chain.values()]))
    return out


def main(argv: list[str] | None = None) -> int:
    ap = make_parser(__name__, __doc__)
    ap.add_argument("runs", help="runs made with --modes or --mode-svd (lines_dump.dat)")
    ap.add_argument("material")
    ap.add_argument(
        "n_modes", type=int, nargs="?", default=6, help="phonon modes per line in the cost proxy (default 6)"
    )
    a = ap.parse_args(argv)
    return emit(lambda: analyze(a.runs, a.material, a.n_modes))


if __name__ == "__main__":
    sys.exit(main())
