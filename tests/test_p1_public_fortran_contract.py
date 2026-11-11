"""Public P1 Fortran module and patched-source falsifiers (no compiler)."""

from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path

from fepdmc_toolkit.p1_da import (
    delayed_acceptance_prob,
    native_re_ratio,
    native_re_target,
    stage2_exact_accept,
)

ROOT = Path(__file__).resolve().parents[1]
MODULE = ROOT / "integration" / "p1" / "linear_da_mod.f90"
RECORD = ROOT / "provenance" / "P1_PUBLIC_PATCH.json"
HIST = ROOT / "src" / "fepdmc_toolkit" / "fortran" / "linear_da_mod.f90"

BANNED = (
    "C2_FORCE_A1",
    "P1_FIXTURE",
    "p1_fixture_on",
    "c1_mark_",
    "c2_dump_swap",
    "c2_live_fp",
    "c0_copy_wq",
    "night1",
    "shadow",
)


def test_adapter_module_is_not_the_benchmarked_module():
    assert MODULE.is_file()
    assert HIST.is_file()
    pub = MODULE.read_bytes()
    hist = HIST.read_bytes()
    assert hashlib.sha256(pub).hexdigest() != hashlib.sha256(hist).hexdigest()
    text = pub.decode()
    for b in BANNED:
        assert b not in text, b
    assert "error stop" in text or "stop '" in text
    assert "DMC_Method" in text
    assert "dmc_band" in text
    assert "zeroTMC" in text
    assert "sample_gt" in text
    assert "Gel(" in text and "Dph(" in text
    assert "cal_ek_int" in text
    # exact ell_R must come from P_accept, not a re-sum of Gel/Dph in stage2
    stage2 = text.split("subroutine linear_da_stage2", 1)[1].split("end subroutine", 1)[0]
    assert "P_accept" in stage2
    assert "log(" in stage2
    assert "ell_hat" in stage2
    assert "Gel(" not in stage2
    assert "Dph(" not in stage2


def test_p1_record_schema_and_not_donor():
    rec = json.loads(RECORD.read_text())
    assert rec["adapter"] == "p1_delayed_acceptance"
    assert rec["architecture"] == "C"
    assert rec["historical_donor_equivalence"] == "not_established"
    assert (
        rec["upstream_commit"]
        == json.loads((ROOT / "provenance" / "UPSTREAM_FEP_DMC.json").read_text())["reference_commit"]
    )
    assert rec["module_sha256"] == hashlib.sha256(MODULE.read_bytes()).hexdigest()
    for name in (
        "PUBLIC_PRISTINE",
        "PUBLIC_C0_APPLIED",
        "PUBLIC_P1_APPLIED",
        "PUBLIC_C0_P1_APPLIED",
    ):
        st = rec["states"][name]
        for key in (
            "perturbo-fep-dmc/pert-src/diagMC_JJ_updates.f90",
            "perturbo-fep-dmc/pert-src/makefile",
            "perturbo-fep-dmc/pert-src/diagMC_JJ.f90",
        ):
            assert len(st[key]) == 64, key


def test_da_reciprocal_identity_with_clipped_score():
    clip = math.log(10.0)
    for r in (0.01, 0.2, 0.9, 1.0, 1.7, 12.0, 40.0):
        ell = math.log(r)
        eh = max(-clip, min(clip, ell))
        rhat = math.exp(eh)
        a_xy = delayed_acceptance_prob(1.0, r, rhat)
        a_yx = delayed_acceptance_prob(r, 1.0, 1.0 / rhat)
        assert abs(a_xy / a_yx - r) <= 1e-12 * max(1.0, r)


def test_dropping_ell_hat_correction_is_a_different_kernel():
    ell_r, ell_hat, ran = 0.3, 0.5, 0.9
    ok = stage2_exact_accept(ran, ell_r, ell_hat)
    wrong = math.log(max(ran, 1.0e-300)) < min(0.0, ell_r)
    assert ok is False
    assert wrong is True


def test_abs_complex_is_not_r_native():
    mx, my = 1.0 + 0.8j, 0.5 - 0.9j
    r = native_re_ratio(my, mx)
    wrong = abs(my) / abs(mx)
    assert abs(wrong - r) > 1e-8
    assert abs(native_re_target(mx)) == abs(mx.real)


def test_zero_new_real_part_rejects():
    assert delayed_acceptance_prob(2.0, 0.0, 1.5) == 0.0


def test_r_native_must_include_p_kchange():
    r_mat = native_re_ratio(0.4 + 0.0j, 0.8 + 0.0j)
    pk = 1.7
    r = r_mat * pk
    assert abs(r_mat - 0.5) <= 1e-12
    assert abs(r - 0.85) <= 1e-12
    assert abs(r - r_mat) > 0.1


def test_odd_clip_preserves_reciprocal_cheap_score():
    clip = math.log(10.0)
    for raw in (-4.0, -2.2, -0.1, 0.0, 0.3, 2.2, 4.0):
        eh = max(-clip, min(clip, raw))
        eh_rev = max(-clip, min(clip, -raw))
        assert abs(eh_rev + eh) <= 1e-15
        assert abs(math.exp(eh) * math.exp(eh_rev) - 1.0) <= 1e-12


def test_nonfinite_exact_ratio_fails_closed():
    assert math.isnan(native_re_ratio(1.0 + 0.0j, complex(float("nan"), 0.0)))
    text = MODULE.read_text()
    stage2 = text.split("subroutine linear_da_stage2", 1)[1].split("end subroutine", 1)[0]
    assert "P_accept is nonfinite" in stage2 or "nonfinite" in stage2
    assert "P_accept <= 0.0_dp" in stage2 or "P_accept <= 0" in stage2
    assert "ell_R = log(P_accept)" in stage2 or "ell_R = log(P_accept" in stage2


def test_public_module_domain_and_rng_contract():
    text = MODULE.read_text()
    assert "DMC_Method must be 0" in text
    assert "dmc_band must be 1" in text
    assert "sample_gt must be false" in text
    assert "zeroTMC must be true" in text
    assert "LINEAR_DA must be off or on" in text
    assert "LINEAR_DA_SCORE must be prop" in text
    assert "9142871" in text
    assert "aux_seed_used == 0_8" in text
    assert "ishft" in text
    assert "diagram%seed" not in text
    assert "RANDOM_NUMBER" not in text
    assert "shadow" not in text.lower()
    assert "C2_FORCE_A1" not in text
    assert "P1_FIXTURE" not in text
    ensure = text.split("subroutine linear_da_ensure", 1)[1].split("end subroutine", 1)[0]
    assert "if (trim(linear_da) == 'on')" in ensure
    assert "DMC_Method" in ensure
    assert "zeroTMC" in ensure
