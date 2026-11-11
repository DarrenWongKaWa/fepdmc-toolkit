"""Independent chord-pairing oracle (D05).

Does not import the recursive generator. Enumeration is leftmost-vertex
perfect matchings of 1..2n. Signs are +1 for bosonic phonons on a single
electron line (no fermion loop).
"""

from __future__ import annotations

from functools import lru_cache

Pairing = tuple[tuple[int, int], ...]


def _matchings(vertices: tuple[int, ...]) -> tuple[Pairing, ...]:
    if not vertices:
        return ((),)
    first = vertices[0]
    out: list[Pairing] = []
    for partner in vertices[1:]:
        rest = tuple(v for v in vertices if v not in (first, partner))
        for tail in _matchings(rest):
            out.append(tuple(sorted(((first, partner),) + tail)))
    return tuple(out)


@lru_cache(maxsize=8)
def all_chord_pairings(n_phonons: int) -> tuple[Pairing, ...]:
    if n_phonons < 0:
        raise ValueError("n_phonons < 0")
    if n_phonons > 6:
        raise ValueError("n_phonons>6 requires an explicit resource gate")
    verts = tuple(range(1, 2 * n_phonons + 1))
    pairings = _matchings(verts)
    if len(pairings) != len(set(pairings)):
        raise RuntimeError("oracle produced duplicate pairings")
    return pairings


def double_factorial_odd(n: int) -> int:
    """(2n-1)!! = number of perfect matchings of 2n points."""
    if n < 0:
        raise ValueError
    acc = 1
    for k in range(n):
        acc *= 2 * k + 1
    return acc


def n_open_sequence(pairing: Pairing) -> tuple[int, ...]:
    n_vertices = 2 * len(pairing)
    seq = []
    for cut in range(1, n_vertices):
        seq.append(sum(1 for a, b in pairing if a <= cut < b))
    return tuple(seq)


def is_proper(pairing: Pairing) -> bool:
    seq = n_open_sequence(pairing)
    return bool(seq) and all(n > 0 for n in seq)


def is_connected_single_electron(pairing: Pairing) -> bool:
    """Single electron line covering 1..2n is always connected as a path."""
    _validate(pairing)
    return True


def is_crossing(pairing: Pairing) -> bool:
    edges = list(pairing)
    for i, (a, b) in enumerate(edges):
        for c, d in edges[i + 1 :]:
            if a < c < b < d or c < a < d < b:
                return True
    return False


def has_nesting(pairing: Pairing) -> bool:
    edges = list(pairing)
    for a, b in edges:
        for c, d in edges:
            if (a, b) == (c, d):
                continue
            if a < c < d < b:
                return True
    return False


def topology_class(pairing: Pairing) -> str:
    _validate(pairing)
    if not is_proper(pairing):
        return "reducible"
    if is_crossing(pairing):
        return "crossed"
    if has_nesting(pairing):
        return "nested"
    return "nested"


def max_open(pairing: Pairing) -> int:
    seq = n_open_sequence(pairing)
    return max(seq) if seq else 0


def classify(pairing: Pairing) -> dict:
    _validate(pairing)
    seq = n_open_sequence(pairing)
    return {
        "n_phonons": len(pairing),
        "pairing": pairing,
        "proper": is_proper(pairing),
        "connected": is_connected_single_electron(pairing),
        "crossing": is_crossing(pairing),
        "nested": has_nesting(pairing),
        "topology_class": topology_class(pairing),
        "n_open": seq,
        "max_open": max(seq) if seq else 0,
        "sign": 1,
        "multiplicity": 1,
    }


def phonon_order(pairing: Pairing) -> Pairing:
    return tuple(sorted(pairing, key=lambda e: (e[0], e[1])))


def _validate(pairing: Pairing) -> None:
    n = len(pairing)
    verts = sorted(v for edge in pairing for v in edge)
    if verts != list(range(1, 2 * n + 1)):
        raise ValueError("pairing must cover 1..2n exactly once")
    if any(a >= b for a, b in pairing):
        raise ValueError("each phonon edge must run earlier -> later")


def pairings_for_class(n_phonons: int, topology_class_name: str) -> tuple[Pairing, ...]:
    all_p = all_chord_pairings(n_phonons)
    if topology_class_name == "nested":
        return tuple(p for p in all_p if is_proper(p) and has_nesting(p) and not is_crossing(p)) or tuple(
            p for p in all_p if is_proper(p) and not is_crossing(p)
        )
    if topology_class_name == "crossed":
        return tuple(p for p in all_p if is_proper(p) and is_crossing(p))
    if topology_class_name == "multi_open":
        return tuple(p for p in all_p if max_open(p) >= 2)
    raise ValueError(f"unknown topology_class {topology_class_name}")
