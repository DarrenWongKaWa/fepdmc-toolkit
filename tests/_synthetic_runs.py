"""Deterministic stand-ins for native FEP-DMC run directories.

The analysis modules read the trace files that the patched ``perturbo.x`` writes.
This module writes files of the same layout with closed-form numbers (no random
generator), so that the analyses can be tested without the native stack. The
numbers have no physical meaning.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

MATERIAL = "lif_hole"
TAU_MAX = 4.0
E_BARE = 9.354873187
N_MEAS = 1200
N_CHAINS = 4


def _wave(n: int, mean: float, amp: float, freq: float, phase: float) -> np.ndarray:
    return mean + amp * np.sin(freq * np.arange(n) + phase)


def _lines(k: int, shift: int) -> list[tuple[int, int]]:
    """Internal phonon lines of sample ``k``: nested, sequential or crossing, with distinct end points."""
    nl = (k + shift) % 5
    if k % 3 == 0:
        pairs = [(j + 1, 2 * nl - j) for j in range(nl)]
    elif k % 3 == 1:
        pairs = [(2 * j + 1, 2 * j + 2) for j in range(nl)]
    else:
        pairs = [(j + 1, j + 1 + nl) for j in range(nl)]
    return sorted(pairs)


def write_chain(root: Path, i: int, dump: str = "M", offset: float = 0.0) -> Path:
    """One chain directory with every trace file; ``dump`` selects the lines_dump.dat layout."""
    c = root / MATERIAL / f"chain{i:02d}"
    c.mkdir(parents=True)
    n = N_MEAS
    idx = np.arange(n)
    ph = 0.7 * i + offset
    sign = np.where((idx * 7 + 3 * i) % 5 == 0, -1.0, 1.0)
    order = 1.0 + ((idx * 3 + i) % 23)
    energy = _wave(n, -0.9 + 0.01 * i + 0.02 * offset, 0.25, 0.37, ph) + 0.1 * np.cos(0.011 * idx + ph)
    num = sign * energy * TAU_MAX
    rb = np.column_stack(
        [
            num,
            sign,
            num * 0.98 + 0.004 * np.cos(0.21 * idx + ph),
            sign * _wave(n, 0.95, 0.04, 0.13, ph),
            (idx + i) % 5,
            1 + (idx + 2 * i) % 5,
            order,
            _wave(n, 0.90, 0.08, 0.29, ph),
            _wave(n, 0.85, 0.12, 0.17, ph + 1.0),
            (idx + i) % 9,
            ((idx + i) % 9) - (idx % 2),
        ]
    )
    np.savetxt(
        c / "rb_trace.dat",
        rb,
        header="raw_num raw_den rb_num rb_den n_elig n_win order rho_single rho_tiled m_tiled m_used",
    )
    (c / "rb_summary.dat").write_text(
        "n_meas 1200\nn_sign_mismatch 0\nn_o_mismatch 0\nmax_o_rel 1.0E-14\nn_zero_fibre 0\nn_tiled_capped 3\n"
        "n_mode_rb_mismatch 0\nn_q_mismatch 0\nmax_q_rel 2.0E-13\n"
    )
    e_chain = num.sum() / sign.sum() / TAU_MAX
    (c / "parameter.dat-1").write_text(f"#Etrue = {e_chain:.12E}\n#other = 1.0\n")

    lines = ((idx + i) % 6).astype(float)
    mode = np.column_stack(
        [
            lines,
            _wave(n, 0.95, 0.04, 0.31, ph),
            _wave(n, 0.50, 0.30, 0.07, ph),
            _wave(n, 0.90, 0.05, 0.23, ph),
            _wave(n, 0.85, 0.05, 0.19, ph),
            _wave(n, 0.50, 0.40, 0.05, ph),
            _wave(n, 0.95, 0.04, 0.41, ph),
            _wave(n, 0.94, 0.05, 0.43, ph + 0.5),
            _wave(n, 0.93, 0.05, 0.47, ph + 1.1),
            np.where(idx % 20 == 0, _wave(n, 0.80, 0.05, 0.3, ph), -1.0),
            _wave(n, 0.70, 0.10, 0.27, ph),
            np.minimum(lines, (idx + i) % 4),
            _wave(n, 0.80, 0.10, 0.33, ph),
            (idx + i) % 3,
            _wave(n, 0.95, 0.04, 0.37, ph),
            _wave(n, 0.94, 0.04, 0.39, ph + 0.3),
            _wave(n, 0.88, 0.05, 0.35, ph + 0.6),
        ]
    )
    np.savetxt(c / "mode_trace.dat", mode, header="mode trace")

    nq = 300
    jq = np.arange(nq)
    q = np.column_stack(
        [
            jq,
            1 + (jq + i) % 7,
            jq % 3,
            1 + (jq * 5 + i) % 14,
            _wave(nq, 0.92, 0.07, 0.29, ph),
            _wave(nq, 0.95, 0.04, 0.31, ph),
        ]
    )
    np.savetxt(c / "q_trace.dat", q, header="q trace")

    mode_rb = np.column_stack([num, sign, num * 0.97, sign * _wave(n, 0.94, 0.03, 0.13, ph), lines])
    np.savetxt(c / "mode_rb_trace.dat", mode_rb, header="raw_num raw_den mrb_num mrb_den n_lines")

    pair = np.where(idx % 3 != 0, _wave(n, 0.85, 0.05, 0.21, ph), -1.0)
    svd = np.column_stack(
        [
            lines,
            _wave(n, 0.95, 0.04, 0.31, ph),
            _wave(n, 0.80, 0.10, 0.27, ph),
            _wave(n, 0.95, 0.04, 0.33, ph),
            _wave(n, 0.80, 0.10, 0.29, ph),
            _wave(n, 2.00, 0.90, 0.11, ph),
            _wave(n, 0.50, 0.40, 0.05, ph),
            _wave(n, 0.50, 0.40, 0.06, ph),
            pair,
            np.abs(pair) * 0.99,
            np.abs(pair) * 0.95,
            np.abs(pair) * 0.94,
        ]
    )
    np.savetxt(c / "mode_svd_trace.dat", svd, header="svd trace")

    out = ["# dump"]
    for k in range(60):
        pairs = _lines(k, i)
        if dump == "M":
            out.append(f"M {k * 10} {len(pairs)} {-1.0 if (k + i) % 5 == 0 else 1.0}")
            out += [f"{a} {b} {0.8 + 0.2 * abs(np.sin(0.9 * k + a + ph)):.10f}" for a, b in pairs]
        else:
            out.append(f"M {k * 10} {len(pairs)} 1.0")
            out.append(f"G {k * 10} {len(pairs)}")
            out += [
                f"{a} {b} {0.8 + 0.2 * abs(np.sin(0.9 * k + a + ph)):.10f} {0.7 + 0.3 * abs(np.cos(0.5 * k + b + ph)):.10f}"
                for a, b in pairs
            ]
    (c / "lines_dump.dat").write_text("\n".join(out) + "\n")

    bchain = np.column_stack(
        [
            num * 0.99,
            sign,
            (idx + i) % 4,
            _wave(n, 0.5, 0.4, 0.05, ph),
            _wave(n, 0.5, 0.4, 0.06, ph),
            _wave(n, 0.5, 0.4, 0.07, ph),
            order,
        ]
    )
    np.savetxt(c / "bchain_trace.dat", bchain, header="f_O f_1 S_eff a b c order")
    (c / "bchain_summary.dat").write_text(
        "n_fold_mismatch 0\nmax_fold_rel 1.2E-13\nn_stage2_accepted 1234\nn_steps 200000\n"
    )
    (c / "stdout.log").write_text(
        "   beta =   4.0000000000E+00\n   Ek-mu =   9354.873187(meV)\n"
        f"            Etrue =    {e_chain:.10E}    0.0E+00\n"
        f"            <g/Z> =     {sign.mean():.10E}    0.0E+00\n"
    )
    (c / "stderr.log").write_text(
        f"\nreal\t1m{2.5 + i + offset:.3f}s\nuser\t0m{50.25 + i + 2 * offset:.3f}s\nsys\t0m0.210s\n"
    )
    (c / "run_meta.json").write_text(
        json.dumps({"bands": [1, 3], "hole": True, "T": 500.0, "seed": 20260925 + 7919 * i})
    )
    return c


def write_runs(root: Path) -> dict[str, Path]:
    """Four run sets (two native-like, one grouped-like, one with the SVD dump) and a tree for collect_sign."""
    sets = {
        "a": (root / "runs_a", "M", 0.0),
        "plain": (root / "runs_plain", "M", 0.1),
        "b": (root / "runs_b", "M", 0.2),
        "svd": (root / "runs_svd", "G", 0.3),
    }
    for path, dump, offset in sets.values():
        for i in range(N_CHAINS):
            write_chain(path, i, dump, offset)
    tree = root / "signtree"
    for group, nmcmc, signs in (
        ("groupA/runs", 200, (0.11, 0.12, 0.10)),
        ("groupB", 200, (0.09, 0.095)),
        ("short/runs", 5, (0.5,)),
    ):
        for i, s in enumerate(signs):
            d = tree / group / f"chain{i:02d}"
            d.mkdir(parents=True)
            (d / "sign.dat-1").write_text(f"{s} 0.001 20000\n")
            (d / "diagMC.in").write_text(f"Nmcmc(1e4) Px Py Pz mu se nq maxOrder\n{nmcmc} 0 0 0 0.0 0 50000 500\n")
            (d / "temper.in").write_text("1 T\n500.00 11.5 1.0E+18\n")
            (d / "pert.in").write_text(
                "&perturbo\n prefix = 'lif'\n calc_mode = 'diagmc-EZ'\n DMC_Method = 0\n band_min = 1\n band_max = 3\n/\n"
            )
    return {name: path for name, (path, _, _) in sets.items()} | {"signtree": tree}


def analysis_arguments(runs: dict[str, Path]) -> dict[str, list[str]]:
    """Command-line arguments of every analysis whose output is compared with the reference file."""
    a, plain, b, svd = (str(runs[k]) for k in ("a", "plain", "b", "svd"))
    eb, tau = str(E_BARE), str(TAU_MAX)
    return {
        "analyze_rb": [a, MATERIAL, tau, eb],
        "analyze_modes": [a, MATERIAL],
        "analyze_q": [a, MATERIAL],
        "analyze_mode_rb": [a, plain, MATERIAL, eb, tau],
        "analyze_svd": [svd, MATERIAL],
        "decompose_native": [a, MATERIAL, tau],
        "estimate_groups": [a, MATERIAL, "6"],
        "compare_bchain": [a, plain, b, MATERIAL, eb, tau],
        "compare_between": [plain, b, MATERIAL, eb, tau],
        "summarize_materials": [a],
    }
