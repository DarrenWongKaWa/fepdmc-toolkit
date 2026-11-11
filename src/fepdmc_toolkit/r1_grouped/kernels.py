"""Markov kernels on the finite grouped toy (docs/r1-grouped-toy.md, section 5).

Two independent code paths share only the configuration helpers:

* ``enumerate_*``: every proposal of a state with its exact probability, used
  to build the sparse transition matrix for exact stationarity checks.
* ``draw_*``: direct random proposals, used by the Monte Carlo chain. Tests
  compare their empirical frequencies with the enumerated probabilities.
"""

from __future__ import annotations

import math
from collections import defaultdict

import numpy as np

from .grouped_toy import (
    NATIVE_MOVES,
    canonical,
    first_eligible_members,
    group_members,
    order,
    partner_map,
    repair,
    window_eligible,
)

# ---------------------------------------------------------------------------
# Native proposals: enumerate
# ---------------------------------------------------------------------------


def _lines(C):
    return list(zip(C[1], C[2]))


def _times_to_lines(C):
    """Lines as ((t_a, t_b), q) in absolute grid times."""
    times = C[0]
    return [((times[a], times[b]), q) for (a, b), q in _lines(C)]


def _from_time_lines(time_lines) -> tuple:
    ts = sorted(t for (a, b), _ in time_lines for t in (a, b))
    pos = {t: i for i, t in enumerate(ts)}
    return canonical(ts, [((pos[a], pos[b]), q) for (a, b), q in time_lines])


def enum_add_remove(model, C):
    n = order(C)
    out = defaultdict(float)
    # add branch (prob 1/2)
    free = [x for x in range(1, model.G + 1) if x not in C[0]]
    if n < model.n_max and len(free) >= 2:
        npairs = len(free) * (len(free) - 1) // 2
        base = _times_to_lines(C)
        for i, a in enumerate(free):
            for b in free[i + 1 :]:
                for q in range(model.L):
                    out[_from_time_lines(base + [((a, b), q)])] += 0.5 / npairs / model.L
    else:
        out[C] += 0.5
    # remove branch (prob 1/2)
    if n >= 1:
        base = _times_to_lines(C)
        for i in range(n):
            out[_from_time_lines(base[:i] + base[i + 1 :])] += 0.5 / n
    else:
        out[C] += 0.5
    return out


def swap_result(C, i):
    times, pairing, q = C
    pm = partner_map(pairing)
    if pm[i] == i + 1:
        return C
    new = []
    for (a, b), qq in zip(pairing, q):
        e = [a, b]
        for k in range(2):
            if e[k] == i:
                e[k] = i + 1
            elif e[k] == i + 1:
                e[k] = i
        new.append(((e[0], e[1]), qq))
    return canonical(times, new)


def enum_swap(model, C):
    n2 = len(C[0])
    out = defaultdict(float)
    if n2 < 4:
        out[C] += 1.0
        return out
    for i in range(n2 - 1):
        out[swap_result(C, i)] += 1.0 / (n2 - 1)
    return out


def shift_candidates(model, C, j):
    times = C[0]
    lo = times[j - 1] if j > 0 else 0
    hi = times[j + 1] if j + 1 < len(times) else model.G + 1
    return [x for x in range(lo + 1, hi) if x != times[j]]


def enum_shift(model, C):
    n2 = len(C[0])
    out = defaultdict(float)
    if n2 == 0:
        out[C] += 1.0
        return out
    for j in range(n2):
        cands = shift_candidates(model, C, j)
        if not cands:
            out[C] += 1.0 / n2
            continue
        for x in cands:
            times = list(C[0])
            times[j] = x
            out[(tuple(times), C[1], C[2])] += 1.0 / n2 / len(cands)
    return out


def enum_qchange(model, C):
    n = order(C)
    out = defaultdict(float)
    if n == 0 or model.L < 2:
        out[C] += 1.0
        return out
    for i in range(n):
        for qq in range(model.L):
            if qq == C[2][i]:
                continue
            q = list(C[2])
            q[i] = qq
            out[(C[0], C[1], tuple(q))] += 1.0 / n / (model.L - 1)
    return out


ENUMERATORS = {
    "add_remove": enum_add_remove,
    "swap": enum_swap,
    "shift": enum_shift,
    "qchange": enum_qchange,
}

# ---------------------------------------------------------------------------
# Native proposals: direct draws (independent of the enumerators)
# ---------------------------------------------------------------------------


def draw_add_remove(model, C, rng):
    n = order(C)
    if rng.random() < 0.5:
        free = [x for x in range(1, model.G + 1) if x not in C[0]]
        if n >= model.n_max or len(free) < 2:
            return C
        a, b = sorted(rng.choice(free, size=2, replace=False).tolist())
        q = int(rng.integers(model.L))
        return _from_time_lines(_times_to_lines(C) + [((a, b), q)])
    if n == 0:
        return C
    i = int(rng.integers(n))
    base = _times_to_lines(C)
    return _from_time_lines(base[:i] + base[i + 1 :])


def draw_swap(model, C, rng):
    n2 = len(C[0])
    if n2 < 4:
        return C
    return swap_result(C, int(rng.integers(n2 - 1)))


def draw_shift(model, C, rng):
    n2 = len(C[0])
    if n2 == 0:
        return C
    j = int(rng.integers(n2))
    cands = shift_candidates(model, C, j)
    if not cands:
        return C
    times = list(C[0])
    times[j] = int(cands[int(rng.integers(len(cands)))])
    return (tuple(times), C[1], C[2])


def draw_qchange(model, C, rng):
    n = order(C)
    if n == 0 or model.L < 2:
        return C
    i = int(rng.integers(n))
    qq = int(rng.integers(model.L - 1))
    if qq >= C[2][i]:
        qq += 1
    q = list(C[2])
    q[i] = qq
    return (C[0], C[1], tuple(q))


DRAWERS = {
    "add_remove": draw_add_remove,
    "swap": draw_swap,
    "shift": draw_shift,
    "qchange": draw_qchange,
}

# ---------------------------------------------------------------------------
# Group moves
# ---------------------------------------------------------------------------


def enum_heatbath_A(model, C, absw, index, *, eligible_only=False):
    """Window heat-bath on |Re D|. Window v uniform on 0..2n-4 (Lemma 1: n is
    invariant). ``eligible_only`` is NC4: v uniform on eligible windows only."""
    out = defaultdict(float)
    n2 = len(C[0])
    if n2 < 4:
        out[C] += 1.0
        return out
    if eligible_only:
        vs = [v for v in range(n2 - 3) if window_eligible(C, v)]
        if not vs:
            out[C] += 1.0
            return out
    else:
        vs = list(range(n2 - 3))
    for v in vs:
        if not window_eligible(C, v):
            out[C] += 1.0 / len(vs)
            continue
        members = [repair(C, v, r) for r in range(3)]
        ws = np.array([absw[index[M]] for M in members])
        tot = ws.sum()
        for M, w in zip(members, ws):
            out[M] += (w / tot) / len(vs)
    return out


def draw_heatbath_A(model, C, absw, index, rng):
    n2 = len(C[0])
    if n2 < 4:
        return C
    v = int(rng.integers(n2 - 3))
    if not window_eligible(C, v):
        return C
    members = [repair(C, v, r) for r in range(3)]
    ws = np.array([absw[index[M]] for M in members])
    return members[int(rng.choice(3, p=ws / ws.sum()))]


def enum_refresh_B(model, C):
    """Uniform refresh inside G(C); Gibbs for any weight constant on groups."""
    members = group_members(C)
    out = defaultdict(float)
    for M in members:
        out[M] += 1.0 / len(members)
    return out


def draw_refresh_B(model, C, rng):
    members = group_members(C)
    return members[int(rng.integers(len(members)))]


def enum_refresh_first(model, C):
    """NC3: uniform refresh inside the first-eligible-window 'group' (not a
    partition, so its weight is not constant on the refreshed set)."""
    members = first_eligible_members(C)
    out = defaultdict(float)
    for M in members:
        out[M] += 1.0 / len(members)
    return out


# ---------------------------------------------------------------------------
# Kernel assembly
# ---------------------------------------------------------------------------


class Kernel:
    """K = (1 - p_group) * sum_m p_m MH_m + p_group * GroupMove.

    ``group`` is one of None, "heatbath_A", "heatbath_A_eligible_only" (NC4),
    "refresh_B". Weights ``w`` are the measure's unnormalized weights.
    """

    def __init__(self, model, states, w, group=None, p_group=0.3):
        self.model = model
        self.states = states
        self.index = {C: i for i, C in enumerate(states)}
        self.w = np.asarray(w, dtype=float)
        self.group = group
        self.p_group = p_group if group else 0.0
        self._qcache: dict = {}

    def _q(self, move, C):
        key = (move, C)
        if key not in self._qcache:
            self._qcache[key] = ENUMERATORS[move](self.model, C)
        return self._qcache[key]

    def _mh_accept(self, move, C, Cp, qf):
        wi = self.w[self.index[C]]
        wj = self.w[self.index[Cp]]
        if wj == 0.0:
            return 0.0
        qr = self._q(move, Cp).get(C, 0.0)
        return min(1.0, (wj * qr) / (wi * qf))

    def row(self, C) -> dict:
        """Exact transition probabilities out of C (requires w(C) > 0)."""
        out = defaultdict(float)
        scale = 1.0 - self.p_group
        for move, pm in NATIVE_MOVES:
            for Cp, qf in self._q(move, C).items():
                if Cp == C:
                    out[C] += scale * pm * qf
                    continue
                a = self._mh_accept(move, C, Cp, qf)
                out[Cp] += scale * pm * qf * a
                out[C] += scale * pm * qf * (1.0 - a)
        if self.group:
            if self.group == "refresh_B":
                g = enum_refresh_B(self.model, C)
            elif self.group == "refresh_first":
                g = enum_refresh_first(self.model, C)
            else:
                g = enum_heatbath_A(
                    self.model, C, self.w, self.index, eligible_only=self.group == "heatbath_A_eligible_only"
                )
            for Cp, p in g.items():
                out[Cp] += self.p_group * p
        return out

    def support(self) -> np.ndarray:
        return np.flatnonzero(self.w > 0)

    def sparse_rows(self):
        """COO triplets over the support (w > 0)."""
        rows, cols, vals = [], [], []
        for i in self.support():
            for Cp, p in self.row(self.states[i]).items():
                if p != 0.0:
                    rows.append(i)
                    cols.append(self.index[Cp])
                    vals.append(p)
        return np.array(rows), np.array(cols), np.array(vals)

    # -- direct Monte Carlo step ---------------------------------------------
    def step(self, C, rng):
        if self.group and rng.random() < self.p_group:
            if self.group == "refresh_B":
                return draw_refresh_B(self.model, C, rng)
            if self.group == "heatbath_A":
                return draw_heatbath_A(self.model, C, self.w, self.index, rng)
            raise ValueError("direct stepping not provided for " + self.group)
        u = rng.random()
        acc = 0.0
        move = NATIVE_MOVES[-1][0]
        for name, pm in NATIVE_MOVES:
            acc += pm
            if u < acc:
                move = name
                break
        Cp = DRAWERS[move](self.model, C, rng)
        if Cp == C:
            return C
        qf = self._q(move, C)[Cp]
        if rng.random() < self._mh_accept(move, C, Cp, qf):
            return Cp
        return C


def stationarity_residual(kernel, pi_full):
    """max |(pi P)_j - pi_j| over all states, and max row-sum error."""
    r, c, v = kernel.sparse_rows()
    piP = np.zeros_like(pi_full)
    np.add.at(piP, c, pi_full[r] * v)
    rowsum = np.zeros_like(pi_full)
    np.add.at(rowsum, r, v)
    sup = kernel.support()
    return float(np.max(np.abs(piP - pi_full))), float(np.max(np.abs(rowsum[sup] - 1.0)))


def detailed_balance_residual(kernel, pi_full):
    r, c, v = kernel.sparse_rows()
    flow = {}
    for i, j, p in zip(r, c, v):
        flow[(i, j)] = flow.get((i, j), 0.0) + pi_full[i] * p
    return max((abs(f - flow.get((j, i), 0.0)) for (i, j), f in flow.items()), default=0.0)


def irreducible_on_support(kernel) -> bool:
    r, c, v = kernel.sparse_rows()
    sup = set(kernel.support().tolist())
    adj = defaultdict(list)
    radj = defaultdict(list)
    for i, j, p in zip(r, c, v):
        if p > 0:
            adj[i].append(j)
            radj[j].append(i)
    start = next(iter(sup))

    def reach(graph):
        seen = {start}
        stack = [start]
        while stack:
            x = stack.pop()
            for y in graph[x]:
                if y not in seen:
                    seen.add(y)
                    stack.append(y)
        return seen

    return reach(adj) >= sup and reach(radj) >= sup


def asymptotic_variance(kernel, pi_full, f_full):
    """sigma^2_as of f under P (exact, Poisson equation), per step."""
    import scipy.sparse as sp
    import scipy.sparse.linalg as spla

    sup = kernel.support()
    pos = {int(s): k for k, s in enumerate(sup)}
    r, c, v = kernel.sparse_rows()
    rr = np.array([pos[int(x)] for x in r])
    cc = np.array([pos[int(x)] for x in c])
    n = len(sup)
    P = sp.csr_matrix((v, (rr, cc)), shape=(n, n))
    pi = pi_full[sup]
    f = f_full[sup]
    fbar = f - float(pi @ f)
    A = (sp.identity(n, format="csr") - P).tolil()
    A[0, :] = 0.0
    A[0, 0] = 1.0
    rhs = fbar.copy()
    rhs[0] = 0.0
    g = spla.spsolve(A.tocsc(), rhs)
    g = g - float(pi @ g)
    return float(2.0 * pi @ (fbar * g) - pi @ (fbar * fbar))


def mc_chain(kernel, C0, nsteps, rng, record):
    """Run the direct chain; ``record`` maps state index -> row of values."""
    idx = kernel.index
    out = np.empty((nsteps, record.shape[1]))
    C = C0
    for s in range(nsteps):
        C = kernel.step(C, rng)
        out[s] = record[idx[C]]
    return out


def batch_means(series, nbatch=50):
    nb = len(series) // nbatch
    b = series[: nb * nbatch].reshape(nbatch, nb, -1).mean(axis=1)
    return b


def ratio_with_error(series_num, series_den, nbatch=50):
    """Ratio of means with jackknife-over-batches standard error."""
    bn = batch_means(series_num[:, None], nbatch)[:, 0]
    bd = batch_means(series_den[:, None], nbatch)[:, 0]
    est = bn.mean() / bd.mean()
    jk = np.array([(bn.sum() - bn[i]) / (bd.sum() - bd[i]) for i in range(nbatch)])
    se = math.sqrt((nbatch - 1) / nbatch * np.sum((jk - jk.mean()) ** 2))
    return float(est), float(se)
