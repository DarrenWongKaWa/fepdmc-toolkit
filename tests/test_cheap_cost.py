"""M7: declared static cost model. No end-to-end speedup claim."""

from __future__ import annotations

import pytest

from fepdmc_toolkit.diagram_compiler import build_scalar_ir, build_twoband_ir, compile_cheap
from fepdmc_toolkit.diagram_compiler.cheap_cost import graph_cost, lower_exact_cost


@pytest.mark.parametrize("n", [1, 2, 3])
def test_twoband_cheap_has_fewer_expensive_primitives(n):
    ir = build_twoband_ir(n)
    exact = lower_exact_cost(ir)
    cheap = compile_cheap(ir).cost
    assert cheap.expensive_primitives < exact.expensive_primitives
    assert cheap.n_vertex == 0
    assert cheap.n_matmul == 0
    assert cheap.n_trace == 0
    assert cheap.n_eigh == 0
    assert exact.n_vertex > 0
    assert exact.n_matmul > 0
    assert exact.n_trace == 1
    assert exact.n_eigh > 0


@pytest.mark.parametrize("n", [1, 2, 3, 4])
def test_scalar_policy_is_exact_so_expensive_count_stays_zero(n):
    ir = build_scalar_ir(n)
    exact = lower_exact_cost(ir)
    cheap = compile_cheap(ir).cost
    assert exact.expensive_primitives == 0
    assert cheap.expensive_primitives == 0
    assert cheap.n_vertex == exact.n_vertex == 0
    assert cheap.n_eigh == 0


def test_cost_model_is_the_declared_node_counts_not_wall_time():
    ir = build_twoband_ir(2)
    cost = graph_cost(compile_cheap(ir).dag, expensive_electron=False)
    assert cost.n_matmul == 0
    assert "wall_time" not in cost.to_dict()
    assert "speedup" not in cost.to_dict()
