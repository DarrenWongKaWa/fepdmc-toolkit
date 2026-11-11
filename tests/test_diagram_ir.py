"""M1/M7: hand-built DiagramIR, canonical builders, serialization round-trip."""

from __future__ import annotations

import json

import pytest

from fepdmc_toolkit.diagram_compiler import (
    MODEL_SCALAR,
    MODEL_TWOBAND,
    Diagram,
    DiagramIR,
    Factor,
    Family,
    IrError,
    MomentumForm,
    Variables,
    build_scalar_diagram,
    build_scalar_ir,
    build_twoband_ir,
    chord_slots,
    diagram_id,
    enumerate_pairings,
    open_chord_slots,
    validate_ir,
)
from fepdmc_toolkit.diagram_compiler.ir import KIND_ELECTRON, KIND_PHONON, KIND_VERTEX


def _scalar_family(n: int, family_id: str = "hand") -> Family:
    return Family(
        family_id=family_id,
        model=MODEL_SCALAR,
        order=n,
        bands=1,
        electron_interface="scalar",
        propagator_convention="physical_unshifted",
        normalization_spec="g_over_sqrt_L",
        series="FULL_ALL_PAIRINGS",
        dispersion="xi_k=2t(1-cos k)",
        phonon_spec="einstein_vacuum_T0",
        rule_version="diagram-compiler-v1",
        global_factors=(Factor("prefactor", normalization="g_over_sqrt_L", exponent="2n"),),
    )


def _variables(n: int) -> Variables:
    return Variables(
        external=("k",),
        chord_slots=tuple(f"q[{i}]" for i in range(n)),
        chord_slot_rule="first_opening_rank",
        time_slots=tuple(f"tau[{i}]" for i in range(2 * n)),
        time_rule="vertex_times_strictly_increasing",
        params=("t", "omega", "g", "L"),
    )


def test_hand_n1_scalar_ir_matches_hand_formula_structure():
    """F_1 = c_1 * G(k-q0, tau1-tau0) * D(0,1): exactly one diagram, two factors."""
    factors = (
        Factor(KIND_ELECTRON, interval=(0, 1), open_chord_slots=(0,), k_form=MomentumForm(1, ((0, -1),))),
        Factor(KIND_PHONON, chord_slot=0, tau_open_idx=0, tau_close_idx=1),
    )
    diagram = Diagram(diagram_id(((1, 2),)), ((1, 2),), 1, 1, factors)
    ir = DiagramIR(_scalar_family(1), _variables(1), (diagram,))
    validate_ir(ir)
    assert len(ir.diagrams) == 1
    e, d = ir.diagrams[0].factors
    assert e.interval == (0, 1) and e.open_chord_slots == (0,)
    assert e.k_form == MomentumForm(1, ((0, -1),))
    assert d.chord_slot == 0 and (d.tau_open_idx, d.tau_close_idx) == (0, 1)


def test_scalar_n2_ir_is_the_three_hand_pairings():
    ir = build_scalar_ir(2)
    validate_ir(ir)
    assert [d.pairing for d in ir.diagrams] == [
        ((1, 2), (3, 4)),
        ((1, 3), (2, 4)),
        ((1, 4), (2, 3)),
    ]
    by_id = {d.diagram_id: d for d in ir.diagrams}
    d = by_id["C_(1,3)_(2,4)"]
    electrons = [f for f in d.factors if f.kind == KIND_ELECTRON]
    assert [(f.interval, f.open_chord_slots) for f in electrons] == [
        ((0, 1), (0,)),
        ((1, 2), (0, 1)),
        ((2, 3), (1,)),
    ]
    assert electrons[1].k_form == MomentumForm(1, ((0, -1), (1, -1)))
    phonons = [f for f in d.factors if f.kind == KIND_PHONON]
    assert [(f.chord_slot, f.tau_open_idx, f.tau_close_idx) for f in phonons] == [(0, 0, 2), (1, 1, 3)]


def test_chord_slot_helpers():
    p = ((1, 3), (2, 4))
    assert chord_slots(p) == {(1, 3): 0, (2, 4): 1}
    assert open_chord_slots(p, 1) == (0,)
    assert open_chord_slots(p, 2) == (0, 1)
    assert open_chord_slots(p, 3) == (1,)
    assert enumerate_pairings(2) == (((1, 2), (3, 4)), ((1, 3), (2, 4)), ((1, 4), (2, 3)))


def test_twoband_n1_ir_factor_order():
    ir = build_twoband_ir(1)
    validate_ir(ir)
    d = ir.diagrams[0]
    kinds = [f.kind for f in d.factors]
    assert kinds == [KIND_VERTEX, KIND_ELECTRON, KIND_VERTEX, KIND_PHONON]
    v1, e, v2, ph = d.factors
    assert (v1.vertex, v1.chord_slot, v1.direction) == (1, 0, "emit")
    assert (v2.vertex, v2.chord_slot, v2.direction) == (2, 0, "absorb")
    assert e.interval == (0, 1) and e.k_form == MomentumForm(1, ((0, -1),))
    assert (ph.chord_slot, ph.tau_open_idx, ph.tau_close_idx) == (0, 0, 1)
    assert not ir.family.global_factors


def test_serialization_round_trip_and_determinism():
    for ir in (build_scalar_ir(2), build_twoband_ir(2)):
        text = ir.to_json()
        assert text == ir.to_json()
        revived = DiagramIR.from_json(text)
        assert revived == ir
        validate_ir(revived)
        payload = json.loads(text)
        assert payload["ir_kind"] == "DiagramIR"
        assert payload["object_type"] == "SHARED_X_GROUP"


def test_serialized_ir_contains_no_sample_numbers():
    def _no_floats(x) -> bool:
        if isinstance(x, float):
            return False
        if isinstance(x, dict):
            return all(_no_floats(v) for v in x.values())
        if isinstance(x, list):
            return all(_no_floats(v) for v in x)
        return True

    assert _no_floats(json.loads(build_scalar_ir(3).to_json()))


def test_from_dict_rejects_foreign_payload():
    with pytest.raises(ValueError):
        DiagramIR.from_dict({"ir_kind": "SomethingElse"})


def _mutated_ir(diagram: Diagram, n: int, *, family=None) -> DiagramIR:
    return DiagramIR(family or _scalar_family(n), _variables(n), (diagram,))


def test_validate_ir_rejects_sign_flip():
    d = build_scalar_diagram(((1, 2),))
    bad = Diagram(d.diagram_id, d.pairing, -1, d.multiplicity, d.factors)
    with pytest.raises(IrError):
        validate_ir(_mutated_ir(bad, 1))


def test_validate_ir_rejects_wrong_k_form_sign():
    d = build_scalar_diagram(((1, 2),))
    factors = tuple(
        Factor(
            f.kind,
            interval=f.interval,
            open_chord_slots=f.open_chord_slots,
            k_form=MomentumForm(1, ((s, 1) for s, _ in f.k_form.q_terms)) if f.k_form else None,
            chord_slot=f.chord_slot,
            tau_open_idx=f.tau_open_idx,
            tau_close_idx=f.tau_close_idx,
        )
        for f in d.factors
    )
    bad = Diagram(d.diagram_id, d.pairing, 1, 1, factors)
    with pytest.raises(IrError):
        validate_ir(_mutated_ir(bad, 1))


def test_validate_ir_rejects_broken_pairing_and_multiplicity():
    with pytest.raises(IrError):
        validate_ir(_mutated_ir(Diagram("x", ((1, 2), (2, 3)), 1, 1, ()), 2))
    d = build_scalar_diagram(((1, 2),))
    with pytest.raises(IrError):
        validate_ir(_mutated_ir(Diagram(d.diagram_id, d.pairing, 1, 2, d.factors), 1))


def test_validate_ir_rejects_duplicate_ids():
    d = build_scalar_diagram(((1, 2),))
    ir = DiagramIR(_scalar_family(1), _variables(1), (d, d))
    with pytest.raises(IrError):
        validate_ir(ir)


def test_validate_ir_rejects_missing_phonon():
    d = build_scalar_diagram(((1, 2),))
    factors = tuple(f for f in d.factors if f.kind != KIND_PHONON)
    bad = Diagram(d.diagram_id, d.pairing, 1, 1, factors)
    with pytest.raises(IrError):
        validate_ir(_mutated_ir(bad, 1))


def test_validate_ir_rejects_wrong_model_interface():
    family = Family(
        family_id="x",
        model=MODEL_TWOBAND,
        order=1,
        bands=1,
        electron_interface="scalar",
        propagator_convention="physical_unshifted",
        normalization_spec="g_over_sqrt_L",
        series="FULL_ALL_PAIRINGS",
        dispersion="xi_k=2t(1-cos k)",
        phonon_spec="einstein_vacuum_T0",
        rule_version="diagram-compiler-v1",
        global_factors=(),
    )
    d = build_scalar_diagram(((1, 2),))
    with pytest.raises(IrError):
        validate_ir(DiagramIR(family, _variables(1), (d,)))


@pytest.mark.parametrize("builder", [build_scalar_ir, build_twoband_ir])
@pytest.mark.parametrize(
    "field,value",
    [
        ("dispersion", "xi=-2t*cos(k)"),
        ("propagator_convention", "unit_G_shifted"),
        ("normalization_spec", "extra_L"),
        ("phonon_spec", "thermal"),
        ("rule_version", "unknown"),
    ],
)
def test_reject_changed_family_semantics(builder, field, value):
    ir = builder(2)
    setattr(ir.family, field, value)
    with pytest.raises(IrError):
        validate_ir(ir)


@pytest.mark.parametrize(
    "field,value",
    [
        ("chord_slot_rule", "closing_rank"),
        ("time_rule", "unordered"),
        ("external", ("p",)),
        ("params", ("g",)),
        ("chord_slots", ("q[1]", "q[0]")),
    ],
)
def test_reject_changed_variable_contract(field, value):
    ir = build_scalar_ir(2)
    setattr(ir.variables, field, value)
    with pytest.raises(IrError):
        validate_ir(ir)


@pytest.mark.parametrize("fault", ["momentum", "interval", "duplicate_phonon", "unknown", "unused_field"])
def test_reject_twoband_semantic_mutations(fault):
    ir = build_twoband_ir(2)
    diagram = ir.diagrams[0]
    electron = next(f for f in diagram.factors if f.kind == KIND_ELECTRON)
    if fault == "momentum":
        electron.k_form = MomentumForm(1, ())
    elif fault == "interval":
        electron.interval = (2, 1)
    elif fault == "duplicate_phonon":
        diagram.factors = diagram.factors[:-1] + (diagram.factors[-2],)
    elif fault == "unknown":
        diagram.factors += (Factor("unknown"),)
    else:
        electron.normalization = "ignored"
    with pytest.raises(IrError):
        validate_ir(ir)


@pytest.mark.parametrize("fault", ["unknown_field", "fractional_order", "bool_sign", "bad_schema", "bad_sign"])
def test_deserialization_fails_closed(fault):
    payload = build_scalar_ir(1).to_dict()
    if fault == "unknown_field":
        payload["diagrams"][0]["factors"][0]["mystery"] = 1
    elif fault == "fractional_order":
        payload["family"]["order"] = 1.9
    elif fault == "bool_sign":
        payload["diagrams"][0]["sign"] = True
    elif fault == "bad_schema":
        payload["schema_version"] = "999"
    else:
        payload["diagrams"][0]["sign"] = -1
    with pytest.raises(ValueError):
        DiagramIR.from_dict(payload)


def test_momentum_keys_are_immutable_and_integral():
    form = MomentumForm(1, ((0, -1),))
    with pytest.raises(AttributeError):
        form.k_coeff = -1
    with pytest.raises(ValueError):
        MomentumForm(1, ((0.9, -1),))


@pytest.mark.parametrize("field", ["interval", "chord_slot", "pairing", "bands"])
def test_validation_rejects_mutated_structural_number_types(field):
    ir = build_scalar_ir(1)
    if field == "interval":
        ir.diagrams[0].factors[0].interval = (0.0, 1.0)
    elif field == "chord_slot":
        ir.diagrams[0].factors[1].chord_slot = False
    elif field == "pairing":
        ir.diagrams[0].pairing = ((1.0, 2.0),)
    else:
        ir.family.bands = True
    with pytest.raises(IrError):
        validate_ir(ir)
