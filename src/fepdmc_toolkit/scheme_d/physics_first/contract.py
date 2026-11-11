"""P-A: declared finite Holstein object. Not a native D(C). Tadpole OFF."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

# Frozen synthetic lattice (not written into the L=2 teacher map as new numbers).
# ξ_k = 2t(1-cos k). Linear Holstein. Phonon vacuum in/out.

OBJECT_ID = "PF_G_FIXED_TIME_ORDER_G2n"
SERIES_FULL = "FULL_ALL_PAIRINGS"
SERIES_PROPER = "PROPER_ONE_ELECTRON_LINE"


@dataclass(frozen=True)
class HolsteinRing:
    L: int
    t: float = 1.0
    omega: float = 0.8
    g: float = 0.45
    tau_end: float = 1.4

    def k_grid(self) -> np.ndarray:
        return 2.0 * np.pi * np.arange(self.L) / self.L

    def xi(self, k: float) -> float:
        return 2.0 * self.t * (1.0 - np.cos(k))

    def volume_vertex(self) -> float:
        return self.g / np.sqrt(self.L)


@dataclass(frozen=True)
class TwoBandHolstein:
    """Hermitian 2-band electron + σ_z Holstein vertex. Physical source, not random."""

    L: int
    t: float = 1.0
    omega: float = 0.8
    g: float = 0.45
    delta: float = 0.35
    gap: float = 0.5

    def H_e(self, k: float) -> np.ndarray:
        xi = 2.0 * self.t * (1.0 - np.cos(k))
        return np.array([[xi, self.delta], [self.delta, xi + self.gap]], dtype=np.complex128)

    def M(self, q: float) -> np.ndarray:
        # M(q)=M(-q)=M† ; Holstein opposite-sign bands.
        return self.g / np.sqrt(self.L) * np.array([[1.0, 0.0], [0.0, -1.0]], dtype=np.complex128)

    def phonon_D(self, dtau: float) -> complex:
        return np.exp(-self.omega * abs(dtau))


def contract_rows() -> list[dict]:
    return [
        {
            "field": "H0_electron",
            "value": "xi_k=2t(1-cos k) on Z_L; one electron",
            "source": "contract.py HolsteinRing.xi",
        },
        {
            "field": "H0_phonon",
            "value": "Einstein Omega; vacuum in and out",
            "source": "contract.py HolsteinRing.omega",
        },
        {
            "field": "Hint",
            "value": "(g/sqrt(L)) sum_{k,q} c^dag_{k-q} c_k (b_q + b^dag_{-q})",
            "source": "linear Holstein; tadpole OFF",
        },
        {
            "field": "object",
            "value": OBJECT_ID,
            "source": "fixed ordered vertex times; not integrated G; not EZ D(C)",
        },
        {
            "field": "measure_q",
            "value": "sum_{q in 2pi n/L} with 1/L in each vertex pair (g/sqrt(L))^{2n}",
            "source": "contract.py volume_vertex",
        },
        {
            "field": "phonon_D",
            "value": "exp(-Omega * (tau_close-tau_open)); T=0 vacuum",
            "source": "contract.py",
        },
        {
            "field": "sign",
            "value": "+1 bosonic; no diagnostic minus on crossing",
            "source": "bosonic phonon lines; a minus sign on crossings would be a different model",
        },
        {
            "field": "tadpole",
            "value": "OFF",
            "source": "tadpole diagrams are not part of the family",
        },
        {
            "field": "not_native_D_of_C",
            "value": "true",
            "source": "the object is not the native FEP-DMC diagram weight D(C)",
        },
    ]
