"""Structural PlanSpec adapter with no dependency on the R1 runtime."""

from __future__ import annotations

from .canonical import build_scalar_ir, build_twoband_ir, validate_ir
from .evaluator import Binding
from .ir import MODEL_SCALAR, MODEL_TWOBAND, DiagramIR

_BUILDERS = {
    MODEL_SCALAR: build_scalar_ir,
    MODEL_TWOBAND: build_twoband_ir,
}

_EXPECTED_NORMALIZATION = {
    MODEL_SCALAR: "g_over_sqrt_L",
    MODEL_TWOBAND: "vertex_contains_g_over_sqrt_L",
}

_EXPECTED_INTERFACE = {MODEL_SCALAR: "scalar", MODEL_TWOBAND: "matrix2"}
_EXPECTED_BANDS = {MODEL_SCALAR: 1, MODEL_TWOBAND: 2}
_EXPECTED_VERTEX_PROVIDER = {
    MODEL_SCALAR: "holstein_scalar_absorbed_in_prefactor",
    MODEL_TWOBAND: "holstein_sigmaz",
}


class ImportRejected(ValueError):
    def __init__(self, code: str, message: str) -> None:
        self.code = code
        super().__init__(f"{code}: {message}")


def import_r1_spec(spec) -> DiagramIR:
    """Compile-once frontend: structure only, no sample numbers."""
    if type(spec.n) is not int or not 1 <= spec.n <= 6:
        raise ImportRejected("IR01", f"order n={spec.n} unsupported")
    expected = {
        "rule_version": "r1-symbolic-v1",
        "expected_material_id": "synthetic",
        "expected_gauge_id": "default",
        "expected_data_version": "analytic_holstein_v1",
    }
    for field, value in expected.items():
        if getattr(spec, field) != value:
            raise ImportRejected("IR13", f"unsupported {field}={getattr(spec, field)!r}")
    if spec.model not in _BUILDERS:
        raise ImportRejected("IR02", f"unsupported model {spec.model!r}")
    if spec.ir_object != "SHARED_X_GROUP" or spec.object_type != "PF_G_FIXED_TIME_ORDER_G2n":
        raise ImportRejected("IR03", f"object {spec.ir_object!r}/{spec.object_type!r} out of Stage 1 scope")
    if spec.propagator_convention != "physical_unshifted":
        raise ImportRejected("IR04", f"unsupported propagator convention {spec.propagator_convention!r}")
    if spec.series != "FULL_ALL_PAIRINGS":
        raise ImportRejected("IR05", f"unsupported series {spec.series!r}")
    if spec.normalization_spec != _EXPECTED_NORMALIZATION[spec.model]:
        raise ImportRejected("IR06", f"unexpected normalization {spec.normalization_spec!r}")
    if spec.electron_interface != _EXPECTED_INTERFACE[spec.model]:
        raise ImportRejected("IR07", f"interface {spec.electron_interface!r} does not match model")
    if spec.bands != _EXPECTED_BANDS[spec.model]:
        raise ImportRejected("IR08", f"bands={spec.bands} does not match model")
    if spec.vertex_provider != _EXPECTED_VERTEX_PROVIDER[spec.model]:
        raise ImportRejected("IR09", f"unexpected vertex provider {spec.vertex_provider!r}")
    if spec.phonon_provider != "einstein_vacuum_T0":
        raise ImportRejected("IR10", f"unsupported phonon provider {spec.phonon_provider!r}")
    if spec.dtype != "complex128":
        raise ImportRejected("IR11", f"unsupported dtype {spec.dtype!r}")
    if spec.boundary_type != "none_current_object":
        raise ImportRejected("IR12", f"unsupported boundary {spec.boundary_type!r}")
    ir = _BUILDERS[spec.model](spec.n)
    validate_ir(ir)
    return ir


def import_r1_binding(binding) -> Binding:
    """Values only. Guards run later in validate_binding (BC* codes)."""
    return Binding(
        k=float(binding.k),
        q=tuple(float(x) for x in binding.q),
        tau=tuple(float(x) for x in binding.tau),
        t=float(binding.t),
        omega=float(binding.omega),
        g=float(binding.g),
        L=binding.L,
        delta=float(binding.delta),
        gap=float(binding.gap),
        finite_model=binding.finite_model,
        continuous_extension=binding.continuous_extension,
        material_id=binding.material_id,
        gauge_id=binding.gauge_id,
        data_version=binding.data_version,
    )
