"""Finite grouped model: two measures on one signed toy space.

Not imported by the C0/P1 adapters and not a material calculation. The
formulae implemented here are the ones written in docs/r1-grouped-toy.md;
section numbers below refer to that file.

Configuration space (section 1). A configuration is

    C = (times, pairing, q)

* ``times``   strictly increasing integers in 1..G (2n vertex times, grid units)
* ``pairing`` tuple of lines (a, b), a < b, over vertex positions 0..2n-1,
              sorted by opening position a (= first-opening rank)
* ``q``       one momentum index in Z_L per line, in the same order.

The external electron runs from time 0 to G+1 with momentum k_ext. Physical
time is ``dtau * grid``.

Weight (section 1). D(C) = tr[U_2n V_2n-1 ... V_0 U_0] * prod_lines Dph with
U_j = exp(-H(k_j) dtau * (t_{j+1} - t_j)), k_j = k_ext - sum(q of lines open
across segment j) mod L, and vertex matrices V from :class:`ToyModel`. With
``eta == 0`` the vertex is (g/sqrt L) sigma_z, H and Dph are the Diagram
Compiler two-band family, and a four-vertex window is exactly its n=2 object.
"""

from __future__ import annotations

import itertools
import math
from dataclasses import dataclass
from functools import cached_property

import numpy as np

SIGMA_Z = np.array([[1.0, 0.0], [0.0, -1.0]], dtype=np.complex128)
SIGMA_X = np.array([[0.0, 1.0], [1.0, 0.0]], dtype=np.complex128)
SIGMA_Y = np.array([[0.0, -1j], [1j, 0.0]], dtype=np.complex128)

# Local pairings of a four-vertex window in window-relative positions.
# Local q slots follow first-opening rank inside the window (section 2).
LOCAL_PAIRINGS = (
    ("sequential", ((0, 1), (2, 3))),
    ("nested", ((0, 3), (1, 2))),
    ("crossing", ((0, 2), (1, 3))),
)

# Native move-type probabilities (state independent, section 5).
NATIVE_MOVES = (("add_remove", 0.30), ("swap", 0.30), ("shift", 0.25), ("qchange", 0.15))


@dataclass(frozen=True)
class ToyModel:
    G: int = 8
    L: int = 2
    n_max: int = 4
    k_ext: int = 0
    t: float = 1.0
    omega: float = 0.5
    g: float = 1.2
    delta: float = 0.9
    gap: float = 0.5
    dtau: float = 0.25
    eta: float = 0.0
    theta0: float = 0.7
    zeta: float = 1.0

    def __post_init__(self) -> None:
        if self.G < 1 or self.L < 1 or self.n_max < 0 or 2 * self.n_max > self.G:
            raise ValueError("need G >= 2 n_max, L >= 1")

    # -- primitives -------------------------------------------------------
    def xi(self, k: int) -> float:
        return 2.0 * self.t * (1.0 - math.cos(2.0 * math.pi * (k % self.L) / self.L))

    def hamiltonian(self, k: int) -> np.ndarray:
        x = self.xi(k)
        return np.array([[x, self.delta], [self.delta, x + self.gap]], dtype=np.complex128)

    @cached_property
    def _u_table(self) -> dict:
        table = {}
        for k in range(self.L):
            w, v = np.linalg.eigh(self.hamiltonian(k))
            for d in range(0, self.G + 2):
                table[(k, d)] = (v * np.exp(-w * self.dtau * d)) @ v.conj().T
        return table

    def U(self, k: int, d: int) -> np.ndarray:
        return self._u_table[(k % self.L, d)]

    def dph(self, d: int) -> float:
        return math.exp(-self.omega * self.dtau * d)

    def vertex(self, direction: str, q: int, k_in: int) -> np.ndarray:
        """e-ph vertex with incoming electron momentum k_in.

        emit:   k_in -> k_in - q,  M_e(k, q) = (g/sqrt L)(sz + eta(cos phi sx + sin phi sy)),
                phi(k, q) = theta0 + 2 pi k / L + 0.9 q  (a k-dependent gauge-like phase)
        absorb: k_in -> k_in + q,  M_a(k_in, q) = M_e(k_in + q, q)^dagger
        eta == 0 and zeta == 1 is the Diagram Compiler vertex (g/sqrt L) sz.
        """
        base = self.g / math.sqrt(self.L)
        if self.eta == 0.0:
            return base * self.zeta * SIGMA_Z
        if direction == "emit":
            k = k_in
        elif direction == "absorb":
            k = k_in + q
        else:
            raise ValueError(direction)
        phi = self.theta0 + 2.0 * math.pi * (k % self.L) / self.L + 0.9 * q
        emit = base * (self.zeta * SIGMA_Z + self.eta * (math.cos(phi) * SIGMA_X + math.sin(phi) * SIGMA_Y))
        return emit if direction == "emit" else emit.conj().T


# ---------------------------------------------------------------------------
# Configurations
# ---------------------------------------------------------------------------


def canonical(times, lines_with_q) -> tuple:
    """Sort lines by opening position; q follows its line."""
    ordered = sorted(((min(a, b), max(a, b)), q) for (a, b), q in lines_with_q)
    return (tuple(times), tuple(e for e, _ in ordered), tuple(q for _, q in ordered))


def order(C) -> int:
    return len(C[1])


def enumerate_pairings(m: int) -> list[tuple]:
    """All perfect matchings of 0..m-1 as edge tuples sorted by opener."""

    def rec(pts):
        if not pts:
            yield ()
            return
        a = pts[0]
        for i in range(1, len(pts)):
            for tail in rec(pts[1:i] + pts[i + 1 :]):
                yield ((a, pts[i]),) + tail

    return [tuple(sorted(p)) for p in rec(tuple(range(m)))]


def enumerate_states(model: ToyModel) -> list[tuple]:
    states = []
    for n in range(model.n_max + 1):
        pairings = enumerate_pairings(2 * n)
        for times in itertools.combinations(range(1, model.G + 1), 2 * n):
            for p in pairings:
                for q in itertools.product(range(model.L), repeat=n):
                    states.append((times, p, q))
    return states


def partner_map(pairing) -> dict:
    out = {}
    for a, b in pairing:
        out[a] = b
        out[b] = a
    return out


def segment_momenta(model: ToyModel, C) -> list[int]:
    """k_j for segments j = 0..2n (segment j ends at vertex j)."""
    times, pairing, q = C
    n2 = len(times)
    ks = []
    for j in range(n2 + 1):
        open_q = sum(qq for (a, b), qq in zip(pairing, q) if a < j <= b)
        ks.append((model.k_ext - open_q) % model.L)
    return ks


def vertex_roles(C) -> list[tuple[str, int]]:
    """(direction, q) of each vertex position."""
    times, pairing, q = C
    roles = [None] * len(times)
    for (a, b), qq in zip(pairing, q):
        roles[a] = ("emit", qq)
        roles[b] = ("absorb", qq)
    return roles


def diagram_weight(model: ToyModel, C) -> complex:
    """Native per-diagram D(C) (BOUND_C evaluation)."""
    times, pairing, q = C
    bounds = (0,) + tuple(times) + (model.G + 1,)
    ks = segment_momenta(model, C)
    roles = vertex_roles(C)
    mat = model.U(ks[0], bounds[1] - bounds[0])
    for j in range(len(times)):
        mat = model.vertex(*roles[j], ks[j]) @ mat
        mat = model.U(ks[j + 1], bounds[j + 2] - bounds[j + 1]) @ mat
    phon = 1.0
    for a, b in pairing:
        phon *= model.dph(times[b] - times[a])
    return complex(np.trace(mat) * phon)


# ---------------------------------------------------------------------------
# Windows, re-pairing, groups (sections 2-3)
# ---------------------------------------------------------------------------


def window_eligible(C, v: int) -> bool:
    """Window v..v+3 holds exactly two lines, both entirely inside it."""
    times, pairing, _ = C
    if v < 0 or v + 3 >= len(times):
        return False
    w = range(v, v + 4)
    touching = [e for e in pairing if e[0] in w or e[1] in w]
    return len(touching) == 2 and all(e[0] in w and e[1] in w for e in touching)


def local_pairing_index(C, v: int) -> int:
    local = tuple(sorted((a - v, b - v) for a, b in C[1] if v <= a <= v + 3))
    for i, (_, p) in enumerate(LOCAL_PAIRINGS):
        if p == local:
            return i
    raise ValueError("window not eligible")


def repair(C, v: int, r: int) -> tuple:
    """rho_{v,r}: replace the window-v pairing by local pairing r; q unchanged.

    The two local lines hold consecutive first-opening ranks (Lemma 1), so
    keeping the q tuple fixed keeps local slot a -> q_R, slot b -> q_{R+1}.
    """
    if not window_eligible(C, v):
        raise ValueError("window not eligible")
    times, pairing, q = C
    outside = [e for e in pairing if not (v <= e[0] <= v + 3)]
    local = [(a + v, b + v) for a, b in LOCAL_PAIRINGS[r][1]]
    new_pairing = tuple(sorted(outside + local))
    return (times, new_pairing, q)


def tiled_windows(C) -> list[int]:
    """Eligible windows among the fixed tiling v = 0, 4, 8, ... (section 3)."""
    n2 = len(C[0])
    return [v for v in range(0, n2 - 3, 4) if window_eligible(C, v)]


def group_members(C) -> list[tuple]:
    """G(C): all combinations of re-pairings of eligible tiled windows."""
    ws = tiled_windows(C)
    members = [C]
    for v in ws:
        members = [repair(M, v, r) for M in members for r in range(3)]
    return members


def first_eligible_members(C) -> list[tuple]:
    """NC3 grouping: re-pairings of the lowest eligible window at ANY position."""
    n2 = len(C[0])
    for v in range(0, n2 - 3):
        if window_eligible(C, v):
            return [repair(C, v, r) for r in range(3)]
    return [C]


# ---------------------------------------------------------------------------
# Summed window operator and factorized group sum (Lemmas 2 and 4)
# ---------------------------------------------------------------------------


def window_ops(model: ToyModel, C, v: int) -> list[np.ndarray]:
    """W^(r) for r = 0,1,2 including local vertices, 3 internal segments and
    the two local phonon lines. k_local is the momentum entering the window."""
    times, pairing, q = C
    ks = segment_momenta(model, C)
    k_local = ks[v]
    R = sum(1 for a, _ in pairing if a < v)
    qa, qb = q[R], q[R + 1]
    tw = times[v : v + 4]
    ops = []
    for _, local in LOCAL_PAIRINGS:
        slot_q = {local[0]: qa, local[1]: qb}
        roles = [None] * 4
        for (a, b), qq in slot_q.items():
            roles[a] = ("emit", qq)
            roles[b] = ("absorb", qq)
        mat = model.vertex(*roles[0], k_local)
        for j in range(1, 4):
            open_q = sum(qq for (a, b), qq in slot_q.items() if a < j <= b)
            k_seg = (k_local - open_q) % model.L
            mat = model.U(k_seg, tw[j] - tw[j - 1]) @ mat
            mat = model.vertex(*roles[j], k_seg) @ mat
        phon = 1.0
        for a, b in local:
            phon *= model.dph(tw[b] - tw[a])
        ops.append(mat * phon)
    return ops


class CompilerWindowOp:
    """R1 participation: the summed window operator from the Diagram Compiler.

    The compiled two-band n=2 SHARED_X_GROUP evaluator returns
    Op = sum_pairings (prod D) M U M U M U M at (k, q_a, q_b, tau[4]); with
    k = k_local this is exactly sum_r W^(r). Valid only for eta == 0.
    """

    def __init__(self, model: ToyModel) -> None:
        if model.eta != 0.0 or model.zeta != 1.0:
            raise ValueError("compiler family has (g/sqrt L) sigma_z vertices only")
        from fepdmc_toolkit.diagram_compiler import Binding, build_twoband_ir, compile_evaluator
        from fepdmc_toolkit.diagram_compiler.evaluator import torus_point

        self._Binding = Binding
        self._torus = torus_point
        self._eval = compile_evaluator(build_twoband_ir(2))
        self.model = model
        self.calls = 0

    def __call__(self, C, v: int) -> np.ndarray:
        m = self.model
        times, pairing, q = C
        k_local = segment_momenta(m, C)[v]
        R = sum(1 for a, _ in pairing if a < v)
        b = self._Binding(
            k=self._torus(k_local, m.L),
            q=(self._torus(q[R], m.L), self._torus(q[R + 1], m.L)),
            tau=tuple(m.dtau * x for x in times[v : v + 4]),
            t=m.t,
            omega=m.omega,
            g=m.g,
            L=m.L,
            delta=m.delta,
            gap=m.gap,
        )
        self.calls += 1
        return np.asarray(self._eval(b)["Op"])


def group_sum(model: ToyModel, C, op_provider=None) -> complex:
    """F(G(C)) by the factorized chain: each eligible tiled window replaced by
    its summed operator (Lemma 4). Never enumerates the 3^m members."""
    times, pairing, q = C
    ws = set(tiled_windows(C))
    bounds = (0,) + tuple(times) + (model.G + 1,)
    ks = segment_momenta(model, C)
    roles = vertex_roles(C)
    local_lines = {e for v in ws for e in pairing if v <= e[0] <= v + 3}
    mat = model.U(ks[0], bounds[1] - bounds[0])
    j = 0
    n2 = len(times)
    while j < n2:
        if j in ws:
            if op_provider is None:
                op = sum(window_ops(model, C, j))
            else:
                op = op_provider(C, j)
            mat = op @ mat
            j += 4
        else:
            mat = model.vertex(*roles[j], ks[j]) @ mat
            j += 1
        mat = model.U(ks[j], bounds[j + 1] - bounds[j]) @ mat
    phon = 1.0
    for a, b in pairing:
        if (a, b) not in local_lines:
            phon *= model.dph(times[b] - times[a])
    return complex(np.trace(mat) * phon)


# ---------------------------------------------------------------------------
# Measures (sections 4, 6)
# ---------------------------------------------------------------------------


class Measure:
    """Unnormalized sampling weight w(C) >= 0 and estimator factor f_O(C) with
    the design identity  sum_C Re D(C) O(C) = sum_C w(C) * f_O(C)."""

    name = "abstract"

    def __init__(self, model: ToyModel, states: list[tuple]) -> None:
        self.model = model
        self.states = states
        self.index = {C: i for i, C in enumerate(states)}
        self.D = np.array([diagram_weight(model, C) for C in states])

    def weights(self) -> np.ndarray:
        raise NotImplementedError

    def estimator(self, obs: np.ndarray) -> np.ndarray:
        raise NotImplementedError


class MeasureA(Measure):
    """BOUND_C: w = |Re D|, f_O = sgn(Re D) O."""

    name = "A_bound_c"

    def weights(self) -> np.ndarray:
        return np.abs(self.D.real)

    def estimator(self, obs):
        return np.sign(self.D.real) * obs


class _GroupedMeasure(Measure):
    grouping = staticmethod(group_members)
    multiplicity = True
    factorized = True

    def __init__(self, model, states, op_provider=None):
        super().__init__(model, states)
        self.group_of = np.empty(len(states), dtype=np.int64)
        self.groups: list[list[int]] = []
        seen: dict[tuple, int] = {}
        for i, C in enumerate(states):
            members = tuple(sorted(self.index[M] for M in self.grouping(C)))
            if members not in seen:
                seen[members] = len(self.groups)
                self.groups.append(list(members))
            self.group_of[i] = seen[members]
        # Explicit member sums (check only) and the consumer's factorized sums.
        self.F_explicit = np.array([self.D[g].sum() for g in self.groups])
        self.op_provider = op_provider
        if self.factorized:
            self.F = np.array([group_sum(model, states[g[0]], op_provider) for g in self.groups])
        else:
            self.F = self.F_explicit

    def weights(self):
        gid = self.group_of
        sizes = np.array([len(g) for g in self.groups])[gid]
        w = np.abs(self.F.real)[gid]
        return w / sizes if self.multiplicity else w

    def estimator(self, obs):
        # f_O(C) = sum_{C' in G} Re D(C') O(C') / |Re F(G)|  (sign times O-bar)
        num = np.array([np.sum(self.D[g].real * obs[g]) for g in self.groups])
        denom = np.abs(self.F.real)
        with np.errstate(divide="ignore", invalid="ignore"):
            per_group = np.where(denom > 0, num / denom, 0.0)
        return per_group[self.group_of]


class MeasureB(_GroupedMeasure):
    """SHARED_X_GROUP on tiled windows: w = |Re F(G)| / |G|."""

    name = "B_grouped"


class NC1NoMultiplicity(_GroupedMeasure):
    """Negative control: grouped weight without the 1/|G| multiplicity."""

    name = "NC1_no_multiplicity"
    multiplicity = False


class NC3FirstEligible(_GroupedMeasure):
    """Negative control: group = re-pairings of the first eligible window at
    any position. Not a partition (Lemma 3 fails)."""

    name = "NC3_first_eligible"
    grouping = staticmethod(first_eligible_members)

    def __init__(self, model, states, op_provider=None):
        Measure.__init__(self, model, states)
        # Each state carries its own (possibly non-closed) member list.
        self.groups = [sorted(self.index[M] for M in first_eligible_members(C)) for C in states]
        self.group_of = np.arange(len(states))
        self.F_explicit = np.array([self.D[g].sum() for g in self.groups])
        self.F = self.F_explicit
        self.op_provider = op_provider


class NC2ModulusWrongSign(Measure):
    """Negative control: w = |D| (complex modulus), reweighted by sgn(Re D)."""

    name = "NC2_modulus_wrong_sign"

    def weights(self):
        return np.abs(self.D)

    def estimator(self, obs):
        return np.sign(self.D.real) * obs


class NC2ModulusCorrect(Measure):
    """Legal counterpart of NC2: w = |D|, reweighted by Re D / |D|."""

    name = "NC2b_modulus_correct"

    def weights(self):
        return np.abs(self.D)

    def estimator(self, obs):
        mod = np.abs(self.D)
        with np.errstate(divide="ignore", invalid="ignore"):
            return np.where(mod > 0, self.D.real / mod, 0.0) * obs


def partition_violations(states, grouping) -> int:
    """Count states C with some C' in grouping(C) whose grouping differs."""
    bad = 0
    for C in states:
        g = set(grouping(C))
        if any(set(grouping(M)) != g for M in g):
            bad += 1
    return bad


# ---------------------------------------------------------------------------
# Observables
# ---------------------------------------------------------------------------


def crossings(C) -> int:
    p = C[1]
    return sum(1 for (a, b), (c, d) in itertools.combinations(p, 2) if a < c < b < d or c < a < d < b)


def observables(states) -> dict[str, np.ndarray]:
    return {
        "order_n": np.array([order(C) for C in states], dtype=float),
        "crossings": np.array([crossings(C) for C in states], dtype=float),
        "one": np.ones(len(states)),
    }
