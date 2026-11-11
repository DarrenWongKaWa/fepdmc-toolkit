"""F03: independent Fock-space time-ordered product. No IR, no pairing list."""

from __future__ import annotations

from itertools import product

import numpy as np

from .contract import HolsteinRing


def _index_k(k: float, L: int) -> int:
    n = int(np.round(k * L / (2.0 * np.pi))) % L
    return n


class FockHolstein:
    """One electron + Einstein phonons on Z_L. Occupations n_q with sum n <= M.

    M >= n_phonon_lines is exact for a 2n-vertex vacuum-to-vacuum insertion
    because at most n phonons are present.
    """

    def __init__(self, ring: HolsteinRing, cutoff: int):
        if cutoff < 0:
            raise ValueError("cutoff")
        self.ring = ring
        self.L = ring.L
        self.M = cutoff
        occ = []
        for ns in product(range(cutoff + 1), repeat=self.L):
            if sum(ns) <= cutoff:
                occ.append(ns)
        self.occ = tuple(occ)
        self.oindex = {o: i for i, o in enumerate(self.occ)}
        self.n_occ = len(self.occ)
        self.dim = self.L * self.n_occ

    def _si(self, k_i: int, o_i: int) -> int:
        return k_i * self.n_occ + o_i

    def energy(self, k_i: int, occ: tuple[int, ...]) -> float:
        k = 2.0 * np.pi * k_i / self.L
        return self.ring.xi(k) + self.ring.omega * sum(occ)

    def apply_V(self, psi: np.ndarray) -> np.ndarray:
        out = np.zeros_like(psi)
        gv = self.ring.volume_vertex()
        L = self.L
        for k_i in range(L):
            for o_i, occ in enumerate(self.occ):
                amp = psi[self._si(k_i, o_i)]
                if amp == 0:
                    continue
                for q in range(L):
                    k_new = (k_i - q) % L
                    qn = (-q) % L
                    # absorb q: b_q
                    if occ[q] > 0:
                        new = list(occ)
                        new[q] -= 1
                        nt = tuple(new)
                        if nt in self.oindex:
                            out[self._si(k_new, self.oindex[nt])] += gv * np.sqrt(occ[q]) * amp
                    # emit -q: b^dag_{-q}
                    if sum(occ) < self.M:
                        new = list(occ)
                        new[qn] += 1
                        nt = tuple(new)
                        if nt in self.oindex:
                            out[self._si(k_new, self.oindex[nt])] += gv * np.sqrt(new[qn]) * amp
        return out

    def propagate(self, psi: np.ndarray, dtau: float) -> np.ndarray:
        out = np.empty_like(psi)
        for k_i in range(self.L):
            for o_i, occ in enumerate(self.occ):
                e = self.energy(k_i, occ)
                out[self._si(k_i, o_i)] = psi[self._si(k_i, o_i)] * np.exp(-e * dtau)
        return out

    def vacuum_to_vacuum(self, k: float, taus: np.ndarray) -> complex:
        """⟨k,0| V G0(Δτ) V ... V |k,0⟩ with 2n vertices at ordered taus."""
        k_i = _index_k(k, self.L)
        n_v = len(taus)
        psi = np.zeros(self.dim, dtype=np.complex128)
        vac = (0,) * self.L
        psi[self._si(k_i, self.oindex[vac])] = 1.0
        psi = self.apply_V(psi)
        for j in range(n_v - 1):
            psi = self.propagate(psi, float(taus[j + 1] - taus[j]))
            psi = self.apply_V(psi)
        return complex(psi[self._si(k_i, self.oindex[vac])])


def fock_group_weight(n: int, taus: np.ndarray, k: float, ring: HolsteinRing) -> complex:
    """Sum over all phonon contractions at fixed times, all q summed inside V.

    This is NOT F_n(x) at one bound q-vector; it is sum_q F_n. Compare only to
    a torus-summed oracle, never to a single-x group.
    """
    fock = FockHolstein(ring, cutoff=n)
    return fock.vacuum_to_vacuum(k, taus)
