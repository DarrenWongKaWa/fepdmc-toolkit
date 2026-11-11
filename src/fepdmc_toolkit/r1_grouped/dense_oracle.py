"""Independent oracle of the grouped toy. Imports nothing from grouped_toy/kernels.

Written from docs/r1-grouped-toy.md, sections 1-4, with different mechanics:

* configurations are generated as sets of time-stamped lines (not position
  pairings), so vertex positions are derived, never stored;
* electron propagators use a scaling-and-squaring Taylor exponential, not an
  eigendecomposition;
* momentum on each segment is recomputed from line time intervals;
* groups are built as equivalence classes of an explicit key (times, q-rank
  vector, lines outside eligible tiled windows) by brute force, and every
  group sum is the explicit sum of its members' per-diagram weights.
"""

from __future__ import annotations

import itertools
import math

import numpy as np


def expm_taylor(A: np.ndarray, terms: int = 30) -> np.ndarray:
    norm = float(np.max(np.sum(np.abs(A), axis=1)))
    s = max(0, int(math.ceil(math.log2(norm))) + 1) if norm > 0.5 else 0
    B = A / (2**s)
    out = np.eye(A.shape[0], dtype=np.complex128)
    term = np.eye(A.shape[0], dtype=np.complex128)
    for k in range(1, terms):
        term = term @ B / k
        out = out + term
    for _ in range(s):
        out = out @ out
    return out


class Oracle:
    def __init__(self, p: dict) -> None:
        self.p = p
        G, L = p["G"], p["L"]
        self.U = {}
        for k in range(L):
            x = 2 * p["t"] * (1 - math.cos(2 * math.pi * k / L))
            H = np.array([[x, p["delta"]], [p["delta"], x + p["gap"]]], dtype=np.complex128)
            for d in range(G + 2):
                self.U[(k, d)] = expm_taylor(-H * p["dtau"] * d)

    def vertex(self, emit: bool, q: int, k_before: int) -> np.ndarray:
        """Emission matrix is written at the electron momentum BEFORE emission;
        absorption is the adjoint of the emission from the momentum AFTER it."""
        p = self.p
        c = p["g"] / math.sqrt(p["L"])
        if p["eta"] == 0.0:
            return c * p["zeta"] * np.diag([1.0, -1.0]).astype(np.complex128)
        k_emit = k_before if emit else k_before + q
        phi = p["theta0"] + 2 * math.pi * (k_emit % p["L"]) / p["L"] + 0.9 * q
        e = p["eta"] * np.exp(1j * phi)
        z = p["zeta"]
        m = c * np.array([[z, np.conj(e)], [e, -z]], dtype=np.complex128)
        return m if emit else m.conj().T

    # ---- configurations as time-stamped line sets --------------------------
    def configurations(self):
        """Yield frozensets of lines (t_open, t_close, rank_q) where rank_q is
        attached by opening time order."""
        p = self.p
        for n in range(p["n_max"] + 1):
            for ts in itertools.combinations(range(1, p["G"] + 1), 2 * n):
                for match in _matchings(list(ts)):
                    opened = sorted(match)
                    for qs in itertools.product(range(p["L"]), repeat=n):
                        yield tuple((a, b, q) for (a, b), q in zip(opened, qs))

    def weight(self, lines) -> complex:
        p = self.p
        ev = sorted([(a, True, q) for a, b, q in lines] + [(b, False, q) for a, b, q in lines])
        ts = [0] + [e[0] for e in ev] + [p["G"] + 1]
        mat = np.eye(2, dtype=np.complex128)
        for s in range(len(ts) - 1):
            mid = ts[s] + 0.5
            k = p["k_ext"] - sum(q for a, b, q in lines if a < mid < b)
            mat = self.U[(k % p["L"], ts[s + 1] - ts[s])] @ mat
            if s < len(ev):
                _, emit, q = ev[s]
                mat = self.vertex(emit, q, k) @ mat
        ph = math.prod(math.exp(-p["omega"] * p["dtau"] * (b - a)) for a, b, q in lines)
        return complex(np.trace(mat) * ph)

    # ---- grouping -----------------------------------------------------------
    @staticmethod
    def tiled_key(lines):
        """Key identifying the group: times, q in opening order, and all lines
        that are not the two local lines of an eligible tiled window."""
        ts = sorted(t for a, b, q in lines for t in (a, b))
        qs = tuple(q for a, b, q in sorted(lines))
        rest = set(lines)
        tiles = []
        for w in range(0, len(ts) - 3, 4):
            win = set(ts[w : w + 4])
            touching = [l for l in lines if l[0] in win or l[1] in win]
            if len(touching) == 2 and all(l[0] in win and l[1] in win for l in touching):
                tiles.append(w)
                rest -= set(touching)
        rest_no_q = tuple(sorted((a, b) for a, b, q in rest))
        return (tuple(ts), qs, rest_no_q, tuple(tiles))

    def exact(self):
        """Exact sums for the physical signed object and both measures."""
        groups: dict = {}
        tot = {"Z_phys": 0.0, "Z_A": 0.0, "Z_mod": 0.0, "sum_ReD_n": 0.0, "sum_ReD_x": 0.0, "states": 0}
        for lines in self.configurations():
            D = self.weight(lines)
            n = len(lines)
            x = sum(
                1
                for l1, l2 in itertools.combinations(lines, 2)
                if l1[0] < l2[0] < l1[1] < l2[1] or l2[0] < l1[0] < l2[1] < l1[1]
            )
            tot["states"] += 1
            tot["Z_phys"] += D.real
            tot["Z_A"] += abs(D.real)
            tot["Z_mod"] += abs(D)
            tot["sum_ReD_n"] += D.real * n
            tot["sum_ReD_x"] += D.real * x
            k = self.tiled_key(lines)
            g = groups.setdefault(k, [0.0 + 0j, 0])
            g[0] += D
            g[1] += 1
        tot["Z_B"] = sum(abs(F.real) for F, _ in groups.values())
        tot["n_groups"] = len(groups)
        tot["group_sizes"] = sorted({c for _, c in groups.values()})
        tot["mean_n"] = tot["sum_ReD_n"] / tot["Z_phys"]
        tot["mean_crossings"] = tot["sum_ReD_x"] / tot["Z_phys"]
        tot["sign_A"] = tot["Z_phys"] / tot["Z_A"]
        tot["sign_B"] = tot["Z_phys"] / tot["Z_B"]
        return tot, groups


def _matchings(pts):
    if not pts:
        yield []
        return
    a = pts[0]
    for i in range(1, len(pts)):
        for rest in _matchings(pts[1:i] + pts[i + 1 :]):
            yield [(a, pts[i])] + rest
