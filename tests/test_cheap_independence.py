"""M6: cheap evaluation must not call the exact path or historical R1."""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

from fepdmc_toolkit.diagram_compiler import Binding, build_scalar_ir, build_twoband_ir, compile_cheap
from fepdmc_toolkit.diagram_compiler.evaluator import torus_point


def _binding(n: int) -> Binding:
    L = 4
    return Binding(
        k=torus_point(0, L),
        q=tuple(torus_point(i + 1, L) for i in range(n)),
        tau=tuple(0.05 * (i + 1) for i in range(2 * n)),
        t=1.0,
        omega=0.8,
        g=0.5,
        L=L,
    )


def test_cheap_does_not_call_exact_evaluator_or_r1(monkeypatch):
    import fepdmc_toolkit.diagram_compiler.evaluator as ev
    import fepdmc_toolkit.diagram_compiler.graph as graph

    def boom(*_args, **_kwargs):
        raise AssertionError("exact or historical backend invoked")

    monkeypatch.setattr(ev, "evaluate_diagramwise", boom)
    monkeypatch.setattr(ev, "evaluate_dag", boom)
    monkeypatch.setattr(ev, "compile_evaluator", boom)
    monkeypatch.setattr(ev, "electron_value", boom)
    monkeypatch.setattr(ev, "vertex_value", boom)
    monkeypatch.setattr(ev, "factor_value", boom)
    monkeypatch.setattr(graph, "lower_diagramwise", boom)
    monkeypatch.setattr(graph, "lower_grouped", boom)

    for builder, n in ((build_scalar_ir, 2), (build_twoband_ir, 2)):
        cheap = compile_cheap(builder(n))
        out = cheap.evaluate(_binding(n))
        assert out["F_hat"].real > 0
        ell = cheap.transition_score(_binding(n), _binding(n))
        assert ell == 0.0


def test_forcing_cheap_to_use_exact_backend_is_detected(monkeypatch):
    import fepdmc_toolkit.diagram_compiler.cheap_evaluator as cheap_mod

    assert "evaluate_diagramwise" not in cheap_mod.__dict__
    called = {"n": 0}

    def spy(*_args, **_kwargs):
        called["n"] += 1
        raise AssertionError("exact backend invoked")

    monkeypatch.setattr(cheap_mod, "evaluate_diagramwise", spy, raising=False)
    compile_cheap(build_scalar_ir(1)).evaluate(_binding(1))
    assert called["n"] == 0


def test_isolated_process_blocks_scheme_d_and_exact_evaluator(tmp_path):
    ir = build_scalar_ir(2)
    binding = _binding(2)
    ir_path = tmp_path / "ir.json"
    bind_path = tmp_path / "binding.json"
    ir_path.write_text(ir.to_json())
    bind_path.write_text(
        json.dumps(
            {
                "k": binding.k,
                "q": list(binding.q),
                "tau": list(binding.tau),
                "t": binding.t,
                "omega": binding.omega,
                "g": binding.g,
                "L": binding.L,
            }
        )
    )
    program = r"""
import importlib.abc, json, sys
from pathlib import Path
class Block(importlib.abc.MetaPathFinder):
    def find_spec(self, fullname, path=None, target=None):
        if fullname.startswith("fepdmc_toolkit.") and not fullname.startswith("fepdmc_toolkit.diagram_compiler"):
            raise AssertionError("old runtime import forbidden: " + fullname)
sys.meta_path.insert(0, Block())
sys.path.insert(0, sys.argv[1])
from fepdmc_toolkit.diagram_compiler import Binding, DiagramIR, compile_cheap
import fepdmc_toolkit.diagram_compiler.evaluator as ev
import fepdmc_toolkit.diagram_compiler.graph as graph
def boom(*a, **k):
    raise AssertionError("exact backend invoked")
ev.evaluate_diagramwise = boom
ev.evaluate_dag = boom
ev.compile_evaluator = boom
ev.electron_value = boom
ev.vertex_value = boom
graph.lower_diagramwise = boom
graph.lower_grouped = boom
ir = DiagramIR.from_json(Path(sys.argv[2]).read_text())
payload = json.loads(Path(sys.argv[3]).read_text())
out = compile_cheap(ir).evaluate(Binding(**payload))
print(out["F_hat"].real)
"""
    src = Path(__file__).resolve().parents[1] / "src"
    result = subprocess.run(
        [sys.executable, "-I", "-c", program, str(src), str(ir_path), str(bind_path)],
        check=True,
        text=True,
        capture_output=True,
    )
    assert float(result.stdout.strip()) > 0.0
