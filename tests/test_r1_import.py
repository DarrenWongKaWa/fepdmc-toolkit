"""M2: R1 frontend adapter. Spec import is structural and fail-closed."""

from __future__ import annotations

from dataclasses import replace

import numpy as np
import pytest

from fepdmc_toolkit.diagram_compiler import (
    Binding,
    BindingError,
    build_scalar_ir,
    build_twoband_ir,
    enumerate_pairings,
    import_r1_binding,
    import_r1_spec,
    validate_binding,
)
from fepdmc_toolkit.scheme_d.oracle import all_chord_pairings
from fepdmc_toolkit.scheme_d.physics_first.symbolic.bindings import base_binding
from fepdmc_toolkit.scheme_d.physics_first.symbolic.plan_spec import scalar_spec, twoband_spec


@pytest.mark.parametrize("n", [1, 2, 3])
def test_scalar_spec_import_equals_canonical_ir(n):
    assert import_r1_spec(scalar_spec(n)) == build_scalar_ir(n)


@pytest.mark.parametrize("n", [1, 2])
def test_twoband_spec_import_equals_canonical_ir(n):
    assert import_r1_spec(twoband_spec(n)) == build_twoband_ir(n)


@pytest.mark.parametrize("n", [1, 2, 3, 4, 5])
def test_pairing_enumeration_matches_scheme_d_oracle(n):
    assert enumerate_pairings(n) == all_chord_pairings(n)


def _rejected(spec):
    with pytest.raises(Exception) as e:
        import_r1_spec(spec)
    assert hasattr(e.value, "code")


def test_import_rejects_unsupported_spec_fields():
    base = scalar_spec(2)
    _rejected(replace(base, n=0))
    _rejected(replace(base, model="Holstein_unknown"))
    _rejected(replace(base, ir_object="BOUND_C"))
    _rejected(replace(base, object_type="SOMETHING_ELSE"))
    _rejected(replace(base, propagator_convention="unit_G_shifted"))
    _rejected(replace(base, series="PROPER_ONE_ELECTRON_LINE"))
    _rejected(replace(base, normalization_spec="g"))
    _rejected(replace(base, electron_interface="matrix2"))
    _rejected(replace(base, bands=2))
    _rejected(replace(base, phonon_provider="free_spectral"))
    _rejected(replace(base, dtype="float64"))
    _rejected(replace(base, boundary_type="open"))
    _rejected(replace(twoband_spec(1), vertex_provider="holstein_scalar_absorbed_in_prefactor"))


def test_binding_import_preserves_values_exactly():
    spec = scalar_spec(2)
    rt = base_binding(spec, "DEV")
    b = import_r1_binding(rt)
    assert b.k == rt.k
    assert b.q == tuple(rt.q)
    assert b.tau == tuple(rt.tau)
    assert (b.t, b.omega, b.g, b.L) == (rt.t, rt.omega, rt.g, rt.L)
    assert b.delta == rt.delta and b.gap == rt.gap
    assert b.finite_model is True


def _binding(**kw):
    base = dict(
        k=0.0,
        q=(0.0, 2.0 * np.pi / 4),
        tau=(0.1, 0.2, 0.3, 0.4),
        t=1.0,
        omega=0.8,
        g=0.45,
        L=4,
    )
    base.update(kw)
    return Binding(**base)


def test_binding_guards_reject_faults():
    ir = import_r1_spec(scalar_spec(2))
    with pytest.raises(BindingError) as e:
        validate_binding(ir, _binding(tau=(0.3, 0.1, 0.2, 0.4)))
    assert e.value.code == "BC01"
    with pytest.raises(BindingError) as e:
        validate_binding(ir, _binding(q=(0.0,)))
    assert e.value.code == "BC02"
    with pytest.raises(BindingError) as e:
        validate_binding(ir, _binding(g=float("nan")))
    assert e.value.code == "BC03"
    with pytest.raises(BindingError) as e:
        validate_binding(ir, _binding(L=0))
    assert e.value.code == "BC04"
    with pytest.raises(BindingError) as e:
        validate_binding(ir, _binding(q=(0.13, 2.0 * np.pi / 4)))
    assert e.value.code == "BC05"


def test_binding_guard_tau_count():
    ir = import_r1_spec(scalar_spec(2))
    with pytest.raises(BindingError) as e:
        validate_binding(ir, _binding(tau=(0.1, 0.2, 0.3)))
    assert e.value.code == "BC02"


@pytest.mark.parametrize(
    "fields",
    [
        {"rule_version": "unknown"},
        {"expected_material_id": "LiF"},
        {"expected_gauge_id": "foreign"},
        {"expected_data_version": "unknown"},
        {"n": 1.9},
        {"n": True},
        {"n": 7},
    ],
)
def test_import_rejects_foreign_contract(fields):
    from fepdmc_toolkit.diagram_compiler.r1_import import ImportRejected

    with pytest.raises(ImportRejected):
        import_r1_spec(replace(scalar_spec(1), **fields))


@pytest.mark.parametrize(
    "fields",
    [
        {"L": 4.9},
        {"L": True},
        {"material_id": "LiF"},
        {"gauge_id": "foreign"},
        {"data_version": "unknown"},
        {"finite_model": True, "continuous_extension": True},
        {"finite_model": False, "continuous_extension": False},
    ],
)
def test_import_does_not_erase_invalid_binding(fields):
    spec = scalar_spec(2)
    rt = replace(base_binding(spec, "DEV"), **fields)
    with pytest.raises((BindingError, ValueError)):
        validate_binding(import_r1_spec(spec), import_r1_binding(rt))


def test_continuous_extension_is_explicit_and_preserved():
    spec = scalar_spec(2)
    rt = replace(base_binding(spec, "DEV"), k=0.123, finite_model=False, continuous_extension=True)
    binding = import_r1_binding(rt)
    assert binding.continuous_extension is True
    validate_binding(import_r1_spec(spec), binding)
