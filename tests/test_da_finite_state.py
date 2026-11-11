"""M6/M7: finite-state detailed balance and stationarity for Design B."""

from __future__ import annotations

import numpy as np

from fepdmc_toolkit.diagram_compiler.da_kernel import (
    transition_matrix,
)


def _balance_residual(pi, p):
    n = len(pi)
    residual = 0.0
    for i in range(n):
        for j in range(n):
            residual = max(residual, abs(pi[i] * p[i, j] - pi[j] * p[j, i]))
    return residual


def test_symmetric_proposal_detailed_balance_and_stationarity():
    pi = np.array([0.5, 0.3, 0.2])
    cheap = np.array([0.1, 0.55, 0.35])
    q = np.array(
        [
            [0.0, 0.5, 0.5],
            [0.5, 0.0, 0.5],
            [0.5, 0.5, 0.0],
        ]
    )
    p = transition_matrix(pi, cheap, q)
    assert _balance_residual(pi, p) <= 1e-12
    left = pi @ p
    assert np.allclose(left, pi, atol=1e-12, rtol=1e-10)
    eig = np.linalg.eig(p.T)
    idx = int(np.argmin(np.abs(eig.eigenvalues - 1.0)))
    vec = np.real(eig.eigenvectors[:, idx])
    vec = vec / vec.sum()
    assert np.allclose(np.abs(vec), pi, atol=1e-8, rtol=1e-6) or np.allclose(vec, pi, atol=1e-8)


def test_asymmetric_reversible_proposal_detailed_balance():
    pi = np.array([0.5, 0.3, 0.2])
    cheap = np.array([0.2, 0.2, 0.6])
    q = np.array(
        [
            [0.0, 0.8, 0.2],
            [0.3, 0.0, 0.7],
            [0.6, 0.4, 0.0],
        ]
    )
    p = transition_matrix(pi, cheap, q)
    assert _balance_residual(pi, p) <= 1e-12
    assert np.allclose(pi @ p, pi, atol=1e-12, rtol=1e-10)


def test_mutation_remove_stage2_breaks_balance():
    pi = np.array([0.5, 0.3, 0.2])
    cheap = np.array([0.1, 0.55, 0.35])
    q = np.array([[0.0, 0.5, 0.5], [0.5, 0.0, 0.5], [0.5, 0.5, 0.0]])
    p = transition_matrix(pi, cheap, q, stage2=False)
    assert _balance_residual(pi, p) > 1e-8


def test_mutation_drop_ell_hat_correction_breaks_balance():
    pi = np.array([0.5, 0.3, 0.2])
    cheap = np.array([0.1, 0.55, 0.35])
    q = np.array([[0.0, 0.5, 0.5], [0.5, 0.0, 0.5], [0.5, 0.5, 0.0]])
    p = transition_matrix(pi, cheap, q, drop_ell_hat_in_stage2=True)
    assert _balance_residual(pi, p) > 1e-8


def test_mutation_invert_hastings_breaks_balance():
    pi = np.array([0.5, 0.3, 0.2])
    cheap = np.array([0.2, 0.2, 0.6])
    q = np.array([[0.0, 0.8, 0.2], [0.3, 0.0, 0.7], [0.6, 0.4, 0.0]])
    p = transition_matrix(pi, cheap, q, invert_hastings=True)
    assert _balance_residual(pi, p) > 1e-8


def test_mutation_cheap_as_exact_breaks_stationarity_on_true_pi():
    pi = np.array([0.5, 0.3, 0.2])
    cheap = np.array([0.1, 0.55, 0.35])
    q = np.array([[0.0, 0.5, 0.5], [0.5, 0.0, 0.5], [0.5, 0.5, 0.0]])
    p = transition_matrix(cheap, cheap, q)
    assert _balance_residual(pi, p) > 1e-8
    assert not np.allclose(pi @ p, pi, atol=1e-8)
