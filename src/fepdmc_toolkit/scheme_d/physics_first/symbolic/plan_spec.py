"""S01: structure spec vs runtime data. Identity excludes sample numbers."""

from __future__ import annotations

from dataclasses import dataclass

from ..contract import OBJECT_ID

RULE_VERSION = "r1-symbolic-v1"
BACKEND = "fortran_ir"
ABI = "packed_complex128_v1"


@dataclass(frozen=True)
class PlanSpec:
    spec_id: str
    object_type: str
    ir_object: str
    model: str
    n: int
    bands: int
    propagator_convention: str
    boundary_type: str
    dtype: str
    rule_version: str
    backend: str
    abi: str
    electron_interface: str
    phonon_provider: str
    vertex_provider: str
    normalization_spec: str
    series: str
    expected_material_id: str = "synthetic"
    expected_gauge_id: str = "default"
    expected_data_version: str = "analytic_holstein_v1"

    def identity_key(self) -> tuple:
        return (
            self.object_type,
            self.ir_object,
            self.model,
            self.n,
            self.bands,
            self.propagator_convention,
            self.boundary_type,
            self.dtype,
            self.rule_version,
            self.backend,
            self.abi,
            self.electron_interface,
            self.phonon_provider,
            self.vertex_provider,
            self.normalization_spec,
            self.series,
        )

    def shape_key(self) -> tuple:
        return (self.n, self.bands, self.electron_interface, self.abi)


def scalar_spec(n: int) -> PlanSpec:
    if n < 1:
        raise ValueError("n>=1")
    return PlanSpec(
        spec_id=f"HS_n{n}_b1",
        object_type=OBJECT_ID,
        ir_object="SHARED_X_GROUP",
        model="Holstein_scalar_vacuum",
        n=n,
        bands=1,
        propagator_convention="physical_unshifted",
        boundary_type="none_current_object",
        dtype="complex128",
        rule_version=RULE_VERSION,
        backend=BACKEND,
        abi=ABI,
        electron_interface="scalar",
        phonon_provider="einstein_vacuum_T0",
        vertex_provider="holstein_scalar_absorbed_in_prefactor",
        normalization_spec="g_over_sqrt_L",
        series="FULL_ALL_PAIRINGS",
    )


def twoband_spec(n: int) -> PlanSpec:
    if n < 1:
        raise ValueError("n>=1")
    return PlanSpec(
        spec_id=f"HT_n{n}_b2",
        object_type=OBJECT_ID,
        ir_object="SHARED_X_GROUP",
        model="Holstein_twoband_hermitian",
        n=n,
        bands=2,
        propagator_convention="physical_unshifted",
        boundary_type="none_current_object",
        dtype="complex128",
        rule_version=RULE_VERSION,
        backend=BACKEND,
        abi=ABI,
        electron_interface="matrix2",
        phonon_provider="einstein_vacuum_T0",
        vertex_provider="holstein_sigmaz",
        normalization_spec="vertex_contains_g_over_sqrt_L",
        series="FULL_ALL_PAIRINGS",
    )


PREREGISTERED_SCALAR = {f"HS_n{n}_b1": scalar_spec(n) for n in (1, 2, 3)}


COUNTEREXAMPLES = (
    ("CE01", "k,q,tau,t,Omega,g,L in node names or plan_id", "runtime_not_identity"),
    ("CE02", "float-equal q or dtau merged as one leaf", "G10"),
    ("CE03", "g or L baked in CONST prefactor", "G11"),
    ("CE04", "builder(spec, taus, qs, k) hidden sample input", "S03"),
    ("CE05", "unit-G skip on raw physical_unshifted", "G06"),
    ("CE06", "silent sort of unsorted times", "G01"),
    ("CE07", "off-grid k/q used as finite torus/Fock evidence", "G05"),
    ("CE08", "n or bands change reused as same spec", "G07"),
    ("CE09", "F_n(x) graph registered as native D(C)", "S15"),
    ("CE10", "old numeric-leaf R1/CSE used as reuse evidence", "S08"),
)
