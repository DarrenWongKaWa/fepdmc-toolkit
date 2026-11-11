"""fepdmc_toolkit.fepdmc.analysis: every analysis on deterministic stand-ins for run directories.

tests/data/analysis_reference.json holds the outputs of the analyses on the files
written by tests/_synthetic_runs.py. The numbers have no physical meaning; they pin
the arithmetic of the modules.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import _synthetic_runs as syn
import pytest

from fepdmc_toolkit.fepdmc import analysis, cli

REFERENCE = json.loads((Path(__file__).resolve().parent / "data" / "analysis_reference.json").read_text())
NEEDS_SCIPY = {"analyze_mode_rb", "compare_between"}


@pytest.fixture(scope="module")
def runs(tmp_path_factory):
    return syn.write_runs(tmp_path_factory.mktemp("runs"))


def assert_same(got, want, where="result"):
    if isinstance(want, dict):
        assert isinstance(got, dict) and sorted(got) == sorted(want), where
        for key in want:
            assert_same(got[key], want[key], f"{where}.{key}")
    elif isinstance(want, list):
        assert isinstance(got, list) and len(got) == len(want), where
        for i, (g, w) in enumerate(zip(got, want)):
            assert_same(g, w, f"{where}[{i}]")
    elif isinstance(want, float) and not isinstance(want, bool):
        if math.isnan(want):
            assert isinstance(got, float) and math.isnan(got), where
        else:
            assert got == pytest.approx(want, rel=1e-9, abs=1e-12), where
    else:
        assert got == want, where


def test_reference_covers_the_listed_analyses():
    assert set(syn.analysis_arguments(dict.fromkeys(("a", "plain", "b", "svd"), Path(".")))) | {"collect_sign"} == set(
        REFERENCE
    )
    assert set(REFERENCE) | {"summarize_sign"} == set(analysis.ANALYSES)


@pytest.mark.parametrize("name", sorted(set(REFERENCE) - {"collect_sign"}))
def test_analysis_reproduces_the_reference_output(name, runs, capsys):
    if name in NEEDS_SCIPY:
        pytest.importorskip("scipy")
    assert analysis.main([name, *syn.analysis_arguments(runs)[name]]) == 0
    assert_same(json.loads(capsys.readouterr().out), REFERENCE[name])


def test_analyses_return_their_result_as_a_dict(runs):
    from fepdmc_toolkit.fepdmc.analysis import analyze_rb, decompose_native

    out = analyze_rb.analyze(str(runs["a"]), syn.MATERIAL, syn.TAU_MAX, syn.E_BARE)
    assert out["chains"] == syn.N_CHAINS
    assert_same(out, REFERENCE["analyze_rb"])
    assert_same(decompose_native.analyze(str(runs["a"]), syn.MATERIAL, syn.TAU_MAX), REFERENCE["decompose_native"])


def test_collect_and_summarize_sign(runs, tmp_path, capsys):
    table = tmp_path / "sign.csv"
    assert analysis.main(["collect_sign", str(runs["signtree"]), str(table)]) == 0
    assert capsys.readouterr().out.strip() == f"6 runs -> {table}"
    assert table.read_text() == REFERENCE["collect_sign"]

    assert analysis.main(["summarize-sign", str(table)]) == 0  # dashes are accepted in the name
    out = json.loads(capsys.readouterr().out)
    assert sorted(out["groups"]) == ["groupA", "groupB"]  # the five-step run is below --min-nmcmc
    assert out["groups"]["groupA"] == {
        "chains": 3,
        "mean_sign": pytest.approx(0.11),
        "stderr": pytest.approx(0.01 / math.sqrt(3)),
    }
    assert out["pooled"]["chains"] == 5
    assert out["pooled"]["mean_sign"] == pytest.approx(0.103)
    assert out["variance_ceiling_perfect_grouping"] == pytest.approx(1 / 0.103**2)

    assert analysis.main(["summarize_sign", str(table), "--min-nmcmc", "1", "--exclude", "groupB"]) == 0
    assert sorted(json.loads(capsys.readouterr().out)["groups"]) == ["groupA", "short"]
    assert analysis.main(["summarize_sign", str(table), "--min-nmcmc", "1000"]) == 2
    assert "no run with at least" in capsys.readouterr().err


def test_summarize_materials_writes_the_optional_file(runs, tmp_path, capsys):
    out = tmp_path / "materials.json"
    assert analysis.main(["summarize_materials", str(runs["a"]), str(out)]) == 0
    assert json.loads(out.read_text()) == json.loads(capsys.readouterr().out)


@pytest.mark.parametrize(
    "argv",
    [
        ["analyze_rb", "{empty}", syn.MATERIAL],
        ["analyze_modes", "{empty}", syn.MATERIAL],
        ["analyze_q", "{empty}", syn.MATERIAL],
        ["analyze_svd", "{empty}", syn.MATERIAL],
        ["decompose_native", "{empty}", syn.MATERIAL],
        ["estimate_groups", "{empty}", syn.MATERIAL],
        ["compare_bchain", "{empty}", "{empty}", "{empty}", syn.MATERIAL, "0.0"],
        ["summarize_materials", "{empty}"],
        ["collect_sign", "{empty}/missing", "{empty}/out.csv"],
        ["summarize_sign", "{empty}/missing.csv"],
    ],
)
def test_missing_input_is_reported_not_raised(argv, tmp_path, capsys):
    (tmp_path / syn.MATERIAL).mkdir()
    argv = [a.replace("{empty}", str(tmp_path)) for a in argv]
    assert analysis.main(argv) == 2
    captured = capsys.readouterr()
    assert captured.err.startswith("error: ") and captured.out == ""


def test_scipy_analyses_report_a_missing_scipy(runs, monkeypatch, capsys):
    import sys

    for mod in list(sys.modules):
        if mod == "scipy" or mod.startswith("scipy."):
            monkeypatch.delitem(sys.modules, mod)
    monkeypatch.setitem(sys.modules, "scipy", None)
    assert analysis.main(["compare_between", *syn.analysis_arguments(runs)["compare_between"]]) == 2
    assert "needs scipy" in capsys.readouterr().err


def test_analyze_command_lists_and_dispatches(capsys):
    assert cli.main(["analyze", "--list"]) == 0
    listed = [line.split()[0] for line in capsys.readouterr().out.splitlines()]
    assert listed == list(analysis.ANALYSES)
    assert cli.main(["analyze", "no_such_analysis"]) == 2
    assert "unknown analysis" in capsys.readouterr().err
    assert cli.main(["analyze"]) == 2
    with pytest.raises(SystemExit) as stop:
        cli.main(["analyze", "analyze_rb", "--help"])
    assert stop.value.code == 0
    assert "fepdmc-toolkit analyze analyze_rb" in capsys.readouterr().out
