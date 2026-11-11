"""SVD mode-basis diagnostic: independence test and scheme estimates.

mode_svd_trace.dat columns:
  n_lines mean_gamma prod_gamma mean_rho prod_rho mean_rank p1 p2
  joint_gamma prod_gamma_pair joint_rho prod_rho_pair
gamma_l = sum_mu|D_mu| / sum_nu|D_nu| (mu: SVD effective modes of the line's
bilinear), rho_l = |sum_nu D_nu| / sum_nu|D_nu|. On single-line fibres both are
exact Z-ratios under pi_A; the pair columns give the exact two-line values.

Schemes (independent-line estimates, validated by the pair test):
  mu basis on all lines:              Z/Z_A ~ E[prod_l gamma_l]
  laminar groups + mu on the rest:    Z/Z_A ~ E[prod_{S} rho_l * prod_{not S} gamma_l]
  laminar groups only (chain B):      Z/Z_A ~ E[prod_{S} rho_l]
"""

import sys
from pathlib import Path

import numpy as np

from ._common import AnalysisError, emit, make_parser
from .estimate_groups import select


def parse_g(path):
    out, cur = [], None
    for line in path.read_text().splitlines():
        if line.startswith("#"):
            continue
        t = line.split()
        if t[0] in ("M", "G"):
            cur = dict(kind=t[0], lines=[])
            out.append(cur)
        elif cur is not None and cur["kind"] == "G":
            cur["lines"].append((int(t[0]), int(t[1]), float(t[2]), float(t[3])))
    return [c for c in out if c["kind"] == "G"]


def analyze(runs: str, mat: str) -> dict:
    res = {"material": mat, "per_chain": {}}
    acc = {
        k: []
        for k in (
            "sign_A",
            "E_prod_gamma",
            "pair_gamma_joint",
            "pair_gamma_prod",
            "pair_rho_joint",
            "pair_rho_prod",
            "E_comb",
            "E_lam",
            "E_mu",
            "E_all_rho",
        )
    }
    for c in sorted((Path(runs) / mat).glob("chain*")):
        a = np.loadtxt(c / "mode_svd_trace.dat")
        r = np.loadtxt(c / "rb_trace.dat")
        sA = float(r[:, 1].mean())
        ok = a[:, 8] >= 0
        rec = dict(
            sign_A=sA,
            mean_rank=float(a[a[:, 0] > 0, 5].mean()),
            E_gamma_line=float(a[a[:, 0] > 0, 1].mean()),
            E_prod_gamma=float(a[:, 2].mean()),
            pair_gamma_joint=float(a[ok, 8].mean()),
            pair_gamma_prod=float(a[ok, 9].mean()),
            pair_rho_joint=float(a[ok, 10].mean()),
            pair_rho_prod=float(a[ok, 11].mean()),
        )
        rec["pair_gamma_joint_over_prod"] = rec["pair_gamma_joint"] / rec["pair_gamma_prod"]
        rec["pair_rho_joint_over_prod"] = rec["pair_rho_joint"] / rec["pair_rho_prod"]
        samples = parse_g(c / "lines_dump.dat")
        comb, lam, mu, allr = [], [], [], []
        for s in samples:
            L = [(x[0], x[1], x[3]) for x in s["lines"]]  # (la, lb, rho) for the selector
            S = set(select(L, 0, 6)) if L else set()
            pr_lam = np.prod([s["lines"][i][3] for i in S]) if S else 1.0
            pr_rest = np.prod([s["lines"][i][2] for i in range(len(L)) if i not in S]) if L else 1.0
            comb.append(pr_lam * pr_rest)
            lam.append(pr_lam)
            mu.append(np.prod([x[2] for x in s["lines"]]) if L else 1.0)
            allr.append(np.prod([x[3] for x in s["lines"]]) if L else 1.0)
        rec.update(
            E_comb=float(np.mean(comb)),
            E_lam=float(np.mean(lam)),
            E_mu=float(np.mean(mu)),
            E_all_rho=float(np.mean(allr)),
            n_dump=len(samples),
        )
        res["per_chain"][c.name] = rec
        for k in acc:
            acc[k].append(rec[k])
    if not res["per_chain"]:
        raise AnalysisError(f"{Path(runs) / mat}: no chain directory")
    m = {k: float(np.mean(v)) for k, v in acc.items()}
    sA = m["sign_A"]
    res["pooled"] = dict(
        sign_A=sA,
        independence_gamma_pair_joint_over_prod=m["pair_gamma_joint"] / m["pair_gamma_prod"],
        independence_rho_pair_joint_over_prod=m["pair_rho_joint"] / m["pair_rho_prod"],
        sign_mu_all_lines=sA / m["E_mu"],
        sign_laminar_groups=sA / m["E_lam"],
        sign_laminar_plus_mu=sA / m["E_comb"],
        sign_all_modes_summed=sA / m["E_all_rho"],
    )
    return res


def main(argv: list[str] | None = None) -> int:
    ap = make_parser(__name__, __doc__)
    ap.add_argument("runs", help="runs made with --mode-svd (mode_svd_trace.dat, lines_dump.dat, rb_trace.dat)")
    ap.add_argument("material")
    a = ap.parse_args(argv)
    return emit(lambda: analyze(a.runs, a.material))


if __name__ == "__main__":
    sys.exit(main())
