"""F02: direct-index weighted oracle. Does not import recurrences or scheme_d helpers."""

from __future__ import annotations

import numpy as np

from .contract import HolsteinRing, TwoBandHolstein

Pairing = tuple[tuple[int, int], ...]


def _matchings(vertices: tuple[int, ...]) -> tuple[Pairing, ...]:
    if not vertices:
        return ((),)
    first = vertices[0]
    out = []
    for partner in vertices[1:]:
        rest = tuple(v for v in vertices if v not in (first, partner))
        for tail in _matchings(rest):
            out.append(tuple(sorted(((first, partner),) + tail)))
    return tuple(out)


def all_pairings(n: int) -> tuple[Pairing, ...]:
    if n < 0 or n > 5:
        raise ValueError("oracle n in 0..5")
    return _matchings(tuple(range(1, 2 * n + 1)))


def is_proper(pairing: Pairing) -> bool:
    n_v = 2 * len(pairing)
    for cut in range(1, n_v):
        if sum(1 for a, b in pairing if a <= cut < b) == 0:
            return False
    return True


def is_crossing(pairing: Pairing) -> bool:
    edges = list(pairing)
    for i, (a, b) in enumerate(edges):
        for c, d in edges[i + 1 :]:
            if a < c < b < d or c < a < d < b:
                return True
    return False


def _g0_scalar(xi: float, dtau: float) -> complex:
    return complex(np.exp(-xi * dtau))


def pairing_weight_scalar(
    pairing: Pairing,
    taus: np.ndarray,
    q_of: dict[tuple[int, int], float],
    k: float,
    ring: HolsteinRing,
    *,
    drop_volume: bool = False,
    flip_phonon_time: bool = False,
    wrong_k: float | None = None,
) -> complex:
    """Explicit product on one electron line. Independent of IR."""
    n = len(pairing)
    verts = 2 * n
    k_use = k if wrong_k is None else wrong_k
    acc = 1.0 + 0j
    gv = 1.0 if drop_volume else ring.volume_vertex()
    for v in range(1, verts + 1):
        acc *= gv
        if v == verts:
            continue
        open_e = [e for e in pairing if e[0] <= v < e[1]]
        k_eff = k_use - sum(q_of[e] for e in open_e)
        dtau = float(taus[v] - taus[v - 1])
        acc *= _g0_scalar(ring.xi(k_eff), dtau)
    for e in pairing:
        dtau_ph = float(taus[e[1] - 1] - taus[e[0] - 1])
        if flip_phonon_time:
            dtau_ph = float(taus[e[0] - 1] - taus[e[1] - 1])
        acc *= np.exp(-ring.omega * dtau_ph)
    return complex(acc)


def q_map_open_order(pairing: Pairing, qs: np.ndarray) -> dict[tuple[int, int], float]:
    """Bind q_j to chords in start-vertex order (dummy labels of x)."""
    ordered = tuple(sorted(pairing, key=lambda e: (e[0], e[1])))
    if len(qs) != len(ordered):
        raise ValueError("qs length")
    return {e: float(qs[i]) for i, e in enumerate(ordered)}


def group_sum_scalar(
    n: int,
    taus: np.ndarray,
    qs: np.ndarray,
    k: float,
    ring: HolsteinRing,
    *,
    series: str = "FULL_ALL_PAIRINGS",
) -> complex:
    acc = 0.0 + 0j
    for p in all_pairings(n):
        if series == "PROPER_ONE_ELECTRON_LINE" and not is_proper(p):
            continue
        acc += pairing_weight_scalar(p, taus, q_map_open_order(p, qs), k, ring)
    return complex(acc)


def _U_eigh(H: np.ndarray, dtau: float) -> np.ndarray:
    w, v = np.linalg.eigh(H)
    return (v * np.exp(-w * dtau)) @ v.conj().T


def pairing_matrix_twoband(
    pairing: Pairing,
    taus: np.ndarray,
    q_of: dict[tuple[int, int], float],
    k: float,
    model: TwoBandHolstein,
    *,
    right_multiply_G: bool = False,
) -> np.ndarray:
    """Untraced 2x2 ordered product. Do not reverse the oracle to match a DAG."""
    n = len(pairing)
    verts = 2 * n
    state = np.eye(2, dtype=np.complex128)
    for v in range(1, verts + 1):
        Mtot = np.eye(2, dtype=np.complex128)
        for e in pairing:
            if e[0] == v or e[1] == v:
                Mtot = model.M(q_of[e]) @ Mtot
        state = Mtot @ state
        if v == verts:
            break
        open_e = [e for e in pairing if e[0] <= v < e[1]]
        k_eff = k - sum(q_of[e] for e in open_e)
        dtau = float(taus[v] - taus[v - 1])
        U = _U_eigh(model.H_e(k_eff), dtau)
        state = state @ U if right_multiply_G else U @ state
    ph = 1.0 + 0j
    for e in pairing:
        dtau_ph = float(taus[e[1] - 1] - taus[e[0] - 1])
        ph *= model.phonon_D(dtau_ph)
    return np.asarray(state * ph, dtype=np.complex128)


def group_sum_matrix_twoband(
    n: int,
    taus: np.ndarray,
    qs: np.ndarray,
    k: float,
    model: TwoBandHolstein,
) -> np.ndarray:
    acc = np.zeros((2, 2), dtype=np.complex128)
    for p in all_pairings(n):
        acc += pairing_matrix_twoband(p, taus, q_map_open_order(p, qs), k, model)
    return acc


def pairing_weight_twoband(
    pairing: Pairing,
    taus: np.ndarray,
    q_of: dict[tuple[int, int], float],
    k: float,
    model: TwoBandHolstein,
    *,
    right_multiply_G: bool = False,
) -> complex:
    return complex(np.trace(pairing_matrix_twoband(pairing, taus, q_of, k, model, right_multiply_G=right_multiply_G)))


def group_sum_twoband(
    n: int,
    taus: np.ndarray,
    qs: np.ndarray,
    k: float,
    model: TwoBandHolstein,
) -> complex:
    acc = 0.0 + 0j
    for p in all_pairings(n):
        acc += pairing_weight_twoband(p, taus, q_map_open_order(p, qs), k, model)
    return complex(acc)


def n1_hand_formula(taus: np.ndarray, q: float, k: float, ring: HolsteinRing) -> complex:
    """n=1 single pairing (1,2). Hand formula, not a loop over matchings."""
    dtau = float(taus[1] - taus[0])
    gv = ring.volume_vertex() ** 2
    g_el = _g0_scalar(ring.xi(k - q), dtau)
    dph = np.exp(-ring.omega * dtau)
    return complex(gv * g_el * dph)


def n2_hand_terms(taus: np.ndarray, qs: np.ndarray, k: float, ring: HolsteinRing) -> dict[str, complex]:
    """Three explicit n=2 pairings. 0-based times; q slots in first-open order.

    Does not call the recurrence or all_pairings.
    Oracle 1-based pairings: (01)(23)->((1,2),(3,4)); (03)(12)->((1,4),(2,3));
    (02)(13)->((1,3),(2,4)).
    """
    if len(taus) != 4 or len(qs) != 2:
        raise ValueError("n=2 hand terms need 4 times and 2 q slots")
    q0, q1 = float(qs[0]), float(qs[1])
    d0 = float(taus[1] - taus[0])
    d1 = float(taus[2] - taus[1])
    d2 = float(taus[3] - taus[2])
    c2 = ring.volume_vertex() ** 4

    def G(p: float, dtau: float) -> complex:
        return _g0_scalar(ring.xi(p), dtau)

    def D(a: int, b: int) -> complex:
        return complex(np.exp(-ring.omega * (float(taus[b]) - float(taus[a]))))

    return {
        "(01)(23)": complex(c2 * G(k - q0, d0) * G(k, d1) * G(k - q1, d2) * D(0, 1) * D(2, 3)),
        "(03)(12)": complex(c2 * G(k - q0, d0) * G(k - q0 - q1, d1) * G(k - q0, d2) * D(0, 3) * D(1, 2)),
        "(02)(13)": complex(c2 * G(k - q0, d0) * G(k - q0 - q1, d1) * G(k - q1, d2) * D(0, 2) * D(1, 3)),
    }


N2_HAND_TO_ORACLE = {
    "(01)(23)": ((1, 2), (3, 4)),
    "(03)(12)": ((1, 4), (2, 3)),
    "(02)(13)": ((1, 3), (2, 4)),
}
