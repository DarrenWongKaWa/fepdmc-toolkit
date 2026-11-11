"""All checks of the grouped toy: exact identities, oracle agreement, stationarity,
negative controls, exact per-step variance, and Monte Carlo confirmation.

Usage:  python -m fepdmc_toolkit.r1_grouped.run_checks [--out DIR] [--steps N] [--regimes R0,R1,R2]

The acceptance thresholds are fixed here (docs/r1-grouped-toy.md, section 8) and
are not tuned after a run. Needs scipy.
"""

from __future__ import annotations

import argparse
import dataclasses
import json
import time
from pathlib import Path

import numpy as np

from .dense_oracle import Oracle
from .grouped_toy import (
    CompilerWindowOp,
    MeasureA,
    MeasureB,
    NC1NoMultiplicity,
    NC2ModulusCorrect,
    NC2ModulusWrongSign,
    NC3FirstEligible,
    ToyModel,
    enumerate_states,
    first_eligible_members,
    group_members,
    observables,
    partition_violations,
    tiled_windows,
)
from .kernels import (
    Kernel,
    asymptotic_variance,
    detailed_balance_residual,
    irreducible_on_support,
    mc_chain,
    ratio_with_error,
    stationarity_residual,
)

TOL_EXACT = 1e-12  # relative, identities and oracle agreement
TOL_STATIONARY = 1e-11  # max|pi P - pi| / max pi
TOL_NC_FAIL = 1e-6  # a negative control must deviate at least this much
MC_SIGMA = 4.0  # MC estimate within 4 standard errors of exact

REGIMES = {
    # Diagram Compiler family: sigma_z vertices; R1 compiled Op is the window operator.
    "R0_compiler": dict(eta=0.0, zeta=1.0, delta=0.9, dtau=0.25, g=1.2, gap=0.5, omega=0.5, k_ext=0, theta0=0.7),
    # Signed: purely interband k-dependent complex vertex.
    "R1_signed_interband": dict(
        eta=1.0, zeta=0.0, delta=2.0, dtau=0.25, g=1.2, gap=0.5, omega=0.2, k_ext=1, theta0=3.46
    ),
    # Signed: mixed vertex, largest grouped sign gain found in the scan.
    "R2_signed_mixed": dict(eta=1.0, zeta=1.0, delta=2.0, dtau=0.4, g=2.0, gap=-1.0, omega=1.0, k_ext=1, theta0=2.63),
}


def rel(a, b):
    a = np.asarray(a)
    b = np.asarray(b)
    return float(np.max(np.abs(a - b)) / max(float(np.max(np.abs(b))), 1e-300))


class Checks:
    def __init__(self):
        self.rows = []

    def add(self, cid, observed, expected, ok, scope):
        self.rows.append(
            dict(check_id=cid, observed=observed, expected=expected, status="PASS" if ok else "FAIL", scope=scope)
        )
        mark = "PASS" if ok else "FAIL"
        print(f"  [{mark}] {cid}: {observed} (expect {expected})", flush=True)

    @property
    def all_pass(self):
        return all(r["status"] == "PASS" for r in self.rows)


def exact_ratio(w, f_num, f_den):
    return float(np.sum(w * f_num) / np.sum(w * f_den))


def run_regime(name, params, steps, seed):
    print(f"== {name} {params}", flush=True)
    model = ToyModel(G=8, L=2, n_max=4, **params)
    states = enumerate_states(model)
    obs = observables(states)
    chk = Checks()
    res = {"regime": name, "params": dataclasses.asdict(model), "n_states": len(states)}

    # ---------------- independent oracle -----------------------------------
    t0 = time.time()
    oracle = Oracle(dataclasses.asdict(model))
    otot, _ = oracle.exact()
    res["oracle_seconds"] = time.time() - t0
    A = MeasureA(model, states)
    ReD = A.D.real
    maxerr = 0.0
    for C, D in zip(states, A.D):
        lines = tuple((C[0][a], C[0][b], q) for (a, b), q in zip(C[1], C[2]))
        maxerr = max(maxerr, abs(oracle.weight(lines) - D))
    scale = float(np.max(np.abs(A.D)))
    chk.add(
        "oracle_D_vs_oracle",
        f"{maxerr / scale:.2e}",
        f"<= {TOL_EXACT}",
        maxerr / scale <= TOL_EXACT,
        "per-diagram D(C): eigh propagators vs Taylor expm, positions vs time-stamped lines",
    )
    chk.add(
        "oracle_state_count",
        otot["states"],
        len(states),
        otot["states"] == len(states),
        "finite configuration space enumeration",
    )
    Zp = float(ReD.sum())
    chk.add(
        "oracle_Z_phys",
        f"{rel(Zp, otot['Z_phys']):.2e}",
        f"<= {TOL_EXACT}",
        rel(Zp, otot["Z_phys"]) <= TOL_EXACT,
        "signed normalization sum Re D",
    )

    provider = CompilerWindowOp(model) if model.eta == 0.0 and model.zeta == 1.0 else None
    B = MeasureB(model, states, op_provider=provider)
    chk.add(
        "identity_lemma3_tiled_partition",
        partition_violations(states, group_members),
        0,
        partition_violations(states, group_members) == 0,
        "tiled-window groups are equivalence classes",
    )
    viol3 = partition_violations(states, first_eligible_members)
    chk.add(
        "NC3_first_eligible_not_partition",
        viol3,
        "> 0",
        viol3 > 0,
        "first-eligible-window grouping breaks closure (negative control)",
    )
    chk.add(
        "oracle_group_count",
        len(B.groups),
        otot["n_groups"],
        len(B.groups) == otot["n_groups"],
        "oracle key-based classes == orbit construction",
    )
    ferr = rel(B.F, B.F_explicit)
    chk.add(
        "identity_lemma4_factorized_group_sum",
        f"{ferr:.2e}",
        f"<= {TOL_EXACT}",
        ferr <= TOL_EXACT,
        "product of summed window operators == explicit member sum"
        + (" (window Op from compiled Diagram Compiler n=2)" if provider else " (own window Op)"),
    )
    ZB = float(np.abs(B.F.real).sum())
    chk.add(
        "oracle_Z_B",
        f"{rel(ZB, otot['Z_B']):.2e}",
        f"<= {TOL_EXACT}",
        rel(ZB, otot["Z_B"]) <= TOL_EXACT,
        "grouped normalization vs oracle",
    )
    if provider:
        res["compiler_window_calls"] = provider.calls
    ZA = float(np.abs(ReD).sum())
    chk.add(
        "identity_sign_theorem",
        f"Z_B/Z_A={ZB / ZA:.6f}",
        "<= 1",
        ZB <= ZA * (1 + 1e-14),
        "Z_B = sum_G |Re F(G)| <= Z_A = sum_C |Re D(C)|",
    )

    m_of = np.array([len(tiled_windows(C)) for C in states])
    absw = np.abs(ReD)
    res.update(
        {
            "Z_phys": Zp,
            "Z_A": ZA,
            "Z_B": ZB,
            "sign_A": Zp / ZA,
            "sign_B": Zp / ZB,
            "sign_gain_B_over_A": ZA / ZB,
            "frac_negative_states": float(np.mean(ReD < 0)),
            "weight_fraction_m_ge_1": float(absw[m_of >= 1].sum() / ZA),
            "group_size_histogram": {
                int(k): int(v) for k, v in zip(*np.unique([len(g) for g in B.groups], return_counts=True))
            },
            "exact_mean_n": otot["mean_n"],
            "exact_mean_crossings": otot["mean_crossings"],
        }
    )

    # ---------------- estimator identities (exact, no sampling) --------------
    measures = {
        "A": A,
        "B": B,
        "NC1": NC1NoMultiplicity(model, states),
        "NC2": NC2ModulusWrongSign(model, states),
        "NC2b": NC2ModulusCorrect(model, states),
        "NC3": NC3FirstEligible(model, states),
    }
    targets = {"order_n": otot["mean_n"], "crossings": otot["mean_crossings"]}
    est = {}
    for mname, M in measures.items():
        w = M.weights()
        f1 = M.estimator(obs["one"])
        est[mname] = {o: exact_ratio(w, M.estimator(obs[o]), f1) for o in targets}
    res["exact_estimator_values"] = est
    for mname in ("A", "B", "NC2b"):
        for o, target in targets.items():
            e = rel(est[mname][o], target)
            chk.add(
                f"identity_unbiased_{mname}_{o}",
                f"{e:.2e}",
                f"<= {TOL_EXACT}",
                e <= TOL_EXACT,
                "E_w[f_O]/E_w[f_1] equals the signed physical expectation",
            )
    complex_regime = float(np.max(np.abs(A.D.imag))) > 1e-12
    for mname, why in (
        ("NC1", "missing 1/|G| multiplicity"),
        ("NC3", "non-partition grouping"),
        ("NC2", "|D| sampled but sgn(Re D) reweighting"),
    ):
        if mname == "NC2" and not complex_regime:
            chk.add(
                "NC2_not_applicable_real_D", "D real", "skip", True, "|D| == |Re D| when D is real; NC2 cannot bias"
            )
            continue
        e = max(rel(est[mname][o], t) for o, t in targets.items())
        chk.add(
            f"{mname}_biased",
            f"{e:.2e}",
            f">= {TOL_NC_FAIL}",
            e >= TOL_NC_FAIL,
            f"negative control ({why}) must be biased",
        )

    # ---------------- exact kernels ---------------------------------------
    kernels = {
        "A_native": Kernel(model, states, A.weights()),
        "A_heatbath": Kernel(model, states, A.weights(), group="heatbath_A"),
        "B_refresh": Kernel(model, states, B.weights(), group="refresh_B"),
        "NC4_A_eligible_only": Kernel(model, states, A.weights(), group="heatbath_A_eligible_only"),
        "NC3_refresh_first": Kernel(model, states, measures["NC3"].weights(), group="refresh_first"),
        "NC1_refresh": Kernel(model, states, measures["NC1"].weights(), group="refresh_B"),
        "B_native": Kernel(model, states, B.weights()),
    }
    res["kernels"] = {}
    for kname, K in kernels.items():
        pi = K.w / K.w.sum()
        st, rowerr = stationarity_residual(K, pi)
        st_rel = st / float(pi.max())
        info = {"stationarity_rel": st_rel, "rowsum_err": rowerr}
        legal = not kname.startswith(("NC4", "NC3"))
        if legal:
            db = detailed_balance_residual(K, pi) / float(pi.max())
            irr = irreducible_on_support(K)
            info.update(detailed_balance_rel=db, irreducible=irr)
            chk.add(
                f"kernel_stationary_{kname}",
                f"{st_rel:.2e}",
                f"<= {TOL_STATIONARY}",
                st_rel <= TOL_STATIONARY and rowerr <= 1e-12,
                "pi P == pi on the finite space",
            )
            chk.add(
                f"kernel_detailed_balance_{kname}",
                f"{db:.2e}",
                f"<= {TOL_STATIONARY}",
                db <= TOL_STATIONARY,
                "every component is reversible",
            )
            chk.add(f"kernel_irreducible_{kname}", irr, True, irr, "support is one communicating class")
        else:
            chk.add(
                f"{kname}_not_stationary",
                f"{st_rel:.2e}",
                f">= {TOL_NC_FAIL}",
                st_rel >= TOL_NC_FAIL,
                "negative control must break pi P == pi",
            )
        res["kernels"][kname] = info

    # Exact per-step asymptotic variance of the ratio estimators (delta method):
    # sigma2 = Var_pi(h) * tau_int with h = (f_O - R f_1) / E[f_1].
    res["asymptotic_variance"] = {}
    res["variance_decomposition"] = {}
    for o in ("order_n", "crossings"):
        res["asymptotic_variance"][o] = {}
        res["variance_decomposition"][o] = {}
        for kname, M in (("A_native", A), ("A_heatbath", A), ("B_native", B), ("B_refresh", B)):
            K = kernels[kname]
            w = K.w / K.w.sum()
            fO = M.estimator(obs[o])
            f1 = M.estimator(obs["one"])
            R = float(w @ fO / (w @ f1))
            h = (fO - R * f1) / float(w @ f1)
            s2 = asymptotic_variance(K, w, h)
            static = float(w @ (h - w @ h) ** 2)
            res["asymptotic_variance"][o][kname] = s2
            res["variance_decomposition"][o][kname] = {"static": static, "tau_int": s2 / static}
        print(
            f"  [{o}] sigma2/step:",
            {k: f"{v:.4g}" for k, v in res["asymptotic_variance"][o].items()},
            " static:",
            {k: f"{v['static']:.4g}" for k, v in res["variance_decomposition"][o].items()},
            " tau_int:",
            {k: f"{v['tau_int']:.3g}" for k, v in res["variance_decomposition"][o].items()},
            flush=True,
        )

    # ---------------- Monte Carlo confirmation ----------------------------
    rng = np.random.default_rng(seed)
    res["mc"] = {}
    for kname, M in (("A_heatbath", A), ("B_refresh", B)):
        K = kernels[kname]
        rec = np.column_stack([M.estimator(obs["one"]), M.estimator(obs["order_n"]), M.estimator(obs["crossings"])])
        start = states[int(np.argmax(K.w))]
        t0 = time.time()
        series = mc_chain(K, start, steps, rng, rec)
        burn = steps // 10
        series = series[burn:]
        out = {"steps": steps, "burn_in": burn, "seconds": time.time() - t0}
        for col, o in ((1, "order_n"), (2, "crossings")):
            e, se = ratio_with_error(series[:, col], series[:, 0])
            z = abs(e - targets[o]) / se
            out[o] = {"estimate": e, "stderr": se, "exact": targets[o], "z": z}
            chk.add(
                f"kernel_mc_{kname}_{o}",
                f"{e:.5f}+-{se:.5f} (z={z:.2f})",
                f"exact {targets[o]:.5f}, z <= {MC_SIGMA}",
                z <= MC_SIGMA,
                "direct-draw MH chain, 50-batch jackknife",
            )
        res["mc"][kname] = out

    res["checks"] = chk.rows
    res["all_pass"] = chk.all_pass
    return res


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="r1_grouped_checks", help="directory for the JSON records")
    ap.add_argument("--steps", type=int, default=400_000)
    ap.add_argument("--seed", type=int, default=20260925)
    ap.add_argument("--regimes", default=",".join(REGIMES))
    a = ap.parse_args(argv)
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    summary = []
    for i, name in enumerate(a.regimes.split(",")):
        r = run_regime(name, REGIMES[name], a.steps, a.seed + i)
        (out / f"{name}.json").write_text(json.dumps(r, indent=2, sort_keys=True, default=float))
        summary.append(
            {
                k: r[k]
                for k in ("regime", "sign_A", "sign_B", "sign_gain_B_over_A", "weight_fraction_m_ge_1", "all_pass")
            }
            | {f"var_{o}_{k}": v for o, d in r["asymptotic_variance"].items() for k, v in d.items()}
        )
    (out / "summary.json").write_text(json.dumps(summary, indent=2, default=float))
    ok = all(s["all_pass"] for s in summary)
    print("ALL PASS" if ok else "SOME CHECKS FAILED")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
