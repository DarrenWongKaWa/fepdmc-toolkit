"""Grouped-measure toy (fepdmc_toolkit.r1_grouped)."""

from __future__ import annotations

import dataclasses
from collections import Counter

import numpy as np
import pytest

from fepdmc_toolkit.r1_grouped.dense_oracle import Oracle
from fepdmc_toolkit.r1_grouped.grouped_toy import (
    CompilerWindowOp,
    MeasureA,
    MeasureB,
    NC1NoMultiplicity,
    NC3FirstEligible,
    ToyModel,
    enumerate_states,
    first_eligible_members,
    group_members,
    group_sum,
    observables,
    partition_violations,
)
from fepdmc_toolkit.r1_grouped.kernels import DRAWERS, ENUMERATORS, Kernel, stationarity_residual

SIGNED = dict(eta=1.0, zeta=0.0, delta=2.0, dtau=0.25, g=1.2, gap=0.5, omega=0.2, k_ext=1, theta0=3.46)


@pytest.fixture(scope="module")
def signed():
    m = ToyModel(G=8, L=2, n_max=4, **SIGNED)
    S = enumerate_states(m)
    return m, S, MeasureA(m, S), MeasureB(m, S)


@pytest.fixture(scope="module")
def small():
    m = ToyModel(G=6, L=2, n_max=3, **SIGNED)
    return m, enumerate_states(m)


def test_weights_match_independent_oracle(signed):
    m, S, A, _ = signed
    o = Oracle(dataclasses.asdict(m))
    got = np.array([o.weight(tuple((C[0][a], C[0][b], q) for (a, b), q in zip(C[1], C[2]))) for C in S])
    assert np.max(np.abs(got - A.D)) <= 1e-12 * np.max(np.abs(A.D))
    assert np.max(np.abs(A.D.imag)) > 1e-6  # the signed regime is genuinely complex
    assert np.any(A.D.real < 0)


def test_tiled_groups_partition_and_first_eligible_does_not(signed):
    _, S, _, _ = signed
    assert partition_violations(S, group_members) == 0
    assert partition_violations(S, first_eligible_members) > 0


def test_factorized_group_sum_equals_member_sum(signed):
    _, _, _, B = signed
    assert np.max(np.abs(B.F - B.F_explicit)) <= 1e-12 * np.max(np.abs(B.F_explicit))
    assert {len(g) for g in B.groups} == {1, 3, 9}


def test_compiler_window_operator_is_the_group_consumer():
    m = ToyModel(G=8, L=2, n_max=4)
    S = enumerate_states(m)
    cw = CompilerWindowOp(m)
    B = MeasureB(m, S, op_provider=cw)
    assert cw.calls > 0
    assert np.max(np.abs(B.F - B.F_explicit)) <= 1e-12 * np.max(np.abs(B.F_explicit))
    with pytest.raises(ValueError):
        CompilerWindowOp(ToyModel(**SIGNED))


def test_sign_theorem_and_estimator_identities(signed):
    m, S, A, B = signed
    obs = observables(S)
    ReD = A.D.real
    ZA, ZB = np.abs(ReD).sum(), np.abs(B.F.real).sum()
    assert ZB < ZA  # strict: some groups cancel in the signed regime
    for o in ("order_n", "crossings"):
        target = np.sum(ReD * obs[o]) / ReD.sum()
        for M in (A, B):
            w = M.weights()
            val = np.sum(w * M.estimator(obs[o])) / np.sum(w * M.estimator(obs["one"]))
            assert abs(val - target) <= 1e-12 * abs(target)
    for NC in (NC1NoMultiplicity(m, S), NC3FirstEligible(m, S)):
        w = NC.weights()
        target = np.sum(ReD * obs["order_n"]) / ReD.sum()
        val = np.sum(w * NC.estimator(obs["order_n"])) / np.sum(w * NC.estimator(obs["one"]))
        assert abs(val - target) > 1e-6


def test_direct_draws_match_enumerated_proposals(small):
    m, S = small
    rng = np.random.default_rng(3)
    C = next(c for c in S if len(c[1]) == 2)
    for move in ENUMERATORS:
        exact = ENUMERATORS[move](m, C)
        N = 40_000
        cnt = Counter(DRAWERS[move](m, C, rng) for _ in range(N))
        assert set(cnt) <= set(exact)
        for Cp, p in exact.items():
            assert abs(cnt[Cp] / N - p) <= 5 * np.sqrt(p * (1 - p) / N) + 1e-9, move


def test_exact_stationarity_legal_and_negative_kernels(small):
    m, S = small
    A = MeasureA(m, S)
    B = MeasureB(m, S)
    for K in (Kernel(m, S, A.weights(), group="heatbath_A"), Kernel(m, S, B.weights(), group="refresh_B")):
        pi = K.w / K.w.sum()
        res, rowerr = stationarity_residual(K, pi)
        assert res <= 1e-13 and rowerr <= 1e-12
    K = Kernel(m, S, A.weights(), group="heatbath_A_eligible_only")
    pi = K.w / K.w.sum()
    assert stationarity_residual(K, pi)[0] > 1e-6


def test_group_sum_never_needs_enumeration(signed):
    m, S, _, B = signed
    C = next(S[g[0]] for g in B.groups if len(g) == 9)
    members = group_members(C)
    assert len(members) == 9
    explicit = sum(MeasureA(m, members).D)
    assert abs(group_sum(m, C) - explicit) <= 1e-12 * abs(explicit)
