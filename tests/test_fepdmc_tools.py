"""fepdmc_toolkit.fepdmc: pooled ratio estimator, exactness validator, patch registry, run driver and CLI."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pytest

from fepdmc_toolkit.fepdmc import pooled as compare_pooled
from fepdmc_toolkit.fepdmc import validate as validate_exactness


def _chain(root: Path, i: int, etrue: float, sign: float, cpu: float) -> None:
    c = root / "lif_hole" / f"chain{i:02d}"
    c.mkdir(parents=True)
    (c / "stdout.log").write_text(
        f"            Etrue =    {etrue:.10E}    0.0E+00\n            <g/Z> =     {sign:.10E}    0.0E+00\n"
    )
    (c / "stderr.log").write_text(f"real\t1m0.0s\nuser\t0m{cpu:.3f}s\nsys\t0m0.1s\n")


def _snapshot(root: Path) -> dict[str, bytes]:
    return {p.relative_to(root).as_posix(): p.read_bytes() for p in sorted(root.rglob("*")) if p.is_file()}


# ---------------------------------------------------------------- pooled estimator


def test_pooled_is_sign_weighted_not_mean_of_ratios(tmp_path):
    # Two chains: a high-sign chain at E=-1 and a low-sign chain at E=-3.
    _chain(tmp_path, 0, -1.0, 0.20, 10.0)
    _chain(tmp_path, 1, -3.0, 0.02, 10.0)
    tau, e_bare = 2.0, 0.0
    side = compare_pooled.side(str(tmp_path), "lif_hole", tau, e_bare)
    expected = (-1.0 * 0.20 + -3.0 * 0.02) / (0.20 + 0.02)
    assert side["Q_pooled"] == pytest.approx(expected)
    assert side["Q_mean_of_chain_ratios"] == pytest.approx(-2.0)
    assert side["cpu_s"] == pytest.approx(10.0)


def test_pooled_needs_two_finished_chains(tmp_path, capsys):
    with pytest.raises(ValueError, match="0 finished chain"):
        compare_pooled.side(str(tmp_path), "lif_hole", 1.0, 0.0)
    _chain(tmp_path, 0, -1.0, 0.2, 1.0)
    (tmp_path / "lif_hole" / "chain01").mkdir()  # started, nothing written yet
    with pytest.raises(ValueError, match="1 finished chain"):
        compare_pooled.side(str(tmp_path), "lif_hole", 1.0, 0.0)
    assert compare_pooled.main([str(tmp_path), str(tmp_path), "lif_hole", "0.0", "1.0"]) == 2
    assert "needs at least two" in capsys.readouterr().err


def test_pooled_tolerates_a_missing_time_log(tmp_path):
    _chain(tmp_path, 0, -1.0, 0.20, 10.0)
    _chain(tmp_path, 1, -1.2, 0.25, 30.0)
    (tmp_path / "lif_hole" / "chain01" / "stderr.log").unlink()
    assert compare_pooled.side(str(tmp_path), "lif_hole", 1.0, 0.0)["cpu_s"] == pytest.approx(10.0)
    (tmp_path / "lif_hole" / "chain00" / "stderr.log").unlink()
    assert np.isnan(compare_pooled.side(str(tmp_path), "lif_hole", 1.0, 0.0)["cpu_s"])


def test_compare_prints_json(tmp_path, capsys):
    rng = np.random.default_rng(1)
    for name, shift in (("a", 0.0), ("c", 0.001)):
        for i in range(4):
            _chain(tmp_path / name, i, -0.87 + shift + 0.001 * rng.standard_normal(), 0.9, 5.0 + i)
    assert compare_pooled.main([str(tmp_path / "a"), str(tmp_path / "c"), "lif_hole", "0.0", "1.0"]) == 0
    out = json.loads(capsys.readouterr().out)
    assert out["A"]["n_chains"] == out["C"]["n_chains"] == 4
    assert out["comparison"]["cpu_ratio_C_over_A"] == pytest.approx(1.0)


# ---------------------------------------------------------------- exactness check


def test_validator_pooled_matches_compare_pooled(tmp_path):
    rng = np.random.default_rng(0)
    for i in range(6):
        _chain(tmp_path, i, -0.87 + 0.001 * rng.standard_normal(), 0.999, 5.0)
    q, se, n = validate_exactness.pooled(tmp_path, "lif_hole")
    side = compare_pooled.side(str(tmp_path), "lif_hole", 1.0, 0.0)
    assert n == 6
    assert q == pytest.approx(side["Q_pooled"])
    assert se == pytest.approx(side["se_delta"])


def test_validator_needs_two_chains(tmp_path):
    _chain(tmp_path, 0, -0.87, 1.0, 1.0)
    q, se, n = validate_exactness.pooled(tmp_path, "lif_hole")
    assert n == 1 and np.isnan(q) and np.isnan(se)


def test_validate_rejects_unknown_kernels(capsys):
    assert validate_exactness.main(["--kernels", "native,bogus", "--analyze-only"]) == 2
    assert validate_exactness.main(["--kernels", "native", "--analyze-only"]) == 2
    assert "at least two of native, pure, mixed" in capsys.readouterr().err


def test_validate_needs_the_dataset_to_start_runs(monkeypatch, capsys):
    monkeypatch.delenv("FEPDMC_TOOLKIT_DATA", raising=False)
    monkeypatch.setattr(validate_exactness.subprocess, "Popen", lambda *a, **k: pytest.fail("no run may start"))
    assert validate_exactness.main([]) == 2
    assert "FEPDMC_TOOLKIT_DATA" in capsys.readouterr().err


def test_validate_analyze_only_compares_existing_runs(tmp_path, capsys):
    ws = tmp_path / "ws"
    args = ["--ws", str(ws), "--analyze-only", "--tag", "v"]
    assert validate_exactness.main(args) == 1
    assert "fewer than two finished chains for: native, pure, mixed" in capsys.readouterr().err

    rng = np.random.default_rng(3)
    for kernel in ("native", "pure", "mixed"):
        for i in range(6):
            _chain(ws / f"v_{kernel}", i, -0.8765 + 0.0005 * rng.standard_normal(), 0.99, 5.0)
    assert validate_exactness.main(args) == 0
    assert capsys.readouterr().out.strip().endswith("PASS")

    for i in range(6):  # a kernel that samples a different target
        (ws / "v_pure" / "lif_hole" / f"chain{i:02d}" / "stdout.log").write_text(
            f"            Etrue =    {-0.8731 + 0.0005 * rng.standard_normal():.10E}    0.0E+00\n"
            "            <g/Z> =     9.9000000000E-01    0.0E+00\n"
        )
    assert validate_exactness.main(args) == 1
    out = capsys.readouterr().out
    assert "DISAGREE" in out and out.strip().endswith("FAIL")


# ---------------------------------------------------------------- patches and prepare


def test_patch_registry_and_data_files():
    from fepdmc_toolkit.fepdmc import FEP_DMC_PIN, patches
    from fepdmc_toolkit.fepdmc.patches import base, bchain

    assert patches.ORDER == ("base", "bchain", "wqfix", "extfix", "extrmfix")
    assert set(patches.PROFILES["fixes"]) == {"base", "wqfix", "extfix", "extrmfix"}
    for name in ("bchain.f90", "rb_window.f90", "mkl_vsl.f90", "make.sys"):
        assert (base.DATA / name).is_file(), name
    assert bchain.DATA == base.DATA
    assert len(FEP_DMC_PIN) == 40
    with pytest.raises(ValueError):
        patches.apply(Path("."), "not-a-patch")


def test_patch_anchor_must_occur_exactly_once():
    from fepdmc_toolkit.fepdmc.patches._anchor import PatchError, locate, once

    assert once("a b c", "b", "x", "label") == "a x c"
    for text in ("a c", "a b b"):
        with pytest.raises(PatchError, match="label"):
            once(text, "b", "x", "label")
    assert locate("abc", "c") == 2
    with pytest.raises(PatchError):
        locate("abc", "z")
    with pytest.raises(PatchError):
        locate("abc", "a", start=1)


def test_pin_matches_provenance_record():
    from fepdmc_toolkit.fepdmc import FEP_DMC_PIN

    rec = json.loads((Path(__file__).resolve().parents[1] / "provenance" / "UPSTREAM_FEP_DMC.json").read_text())
    assert rec["reference_commit"] == FEP_DMC_PIN


def test_prepare_refuses_tree_without_pert_src(tmp_path):
    from fepdmc_toolkit.fepdmc.prepare import main, prepare_tree

    with pytest.raises(FileNotFoundError):
        prepare_tree(tmp_path)
    assert main([str(tmp_path)]) == 2
    (tmp_path / "pert-src").mkdir()
    with pytest.raises(ValueError, match="unknown profile"):
        prepare_tree(tmp_path, "no-such-profile")


def test_prepare_refuses_an_unknown_tree_and_leaves_it_unchanged(tmp_path, capsys):
    from fepdmc_toolkit.fepdmc import PatchError
    from fepdmc_toolkit.fepdmc.prepare import main, prepare_tree

    src = tmp_path / "pert-src"
    src.mkdir()
    (src / "random_tool.f90").write_text("! not the pinned source\n")
    (tmp_path / "make.sys").write_text("the user's make.sys\n")
    before = _snapshot(tmp_path)
    for dry_run in (True, False):
        with pytest.raises(PatchError, match="patch 'base' does not apply"):
            prepare_tree(tmp_path, "fixes", dry_run=dry_run)
        assert _snapshot(tmp_path) == before
    (src / "random_tool.f90").unlink()
    with pytest.raises(PatchError, match="random_tool.f90 is missing"):
        prepare_tree(tmp_path, "research")
    (src / "random_tool.f90").write_text("! not the pinned source\n")
    assert main([str(tmp_path), "--profile", "fixes"]) == 1
    assert "No file was written" in capsys.readouterr().err
    assert _snapshot(tmp_path) == before


# ---------------------------------------------------------------- run driver


def test_runner_switch_env_maps_flags():
    from fepdmc_toolkit.fepdmc import runner

    a = runner.parser().parse_args(
        ["ez", "lif_hole", "2", "--wqfix", "--extfix", "--extrmfix", "--extar", "0.02", "--data", "/d"]
    )
    env, rb = runner.switch_env(a, "lif_hole")
    assert (env["FUTUREB_WQFIX"], env["FUTUREB_EXTFIX"], env["FUTUREB_EXTRMFIX"]) == ("1", "1", "1")
    assert (env["FUTUREB_EXTAR"], env["FUTUREB_EXTAR_P"]) == ("1", "0.02")
    assert env["FUTUREB_ANY"] == "0" and "FUTUREB_MV_UNTIL" not in env
    assert rb == "0"  # three-band hole: no single-band RB
    assert runner.binary("perturbo-fep-dmc-pkg").endswith("perturbo-fep-dmc-pkg/pert-src/perturbo.x")
    assert runner.binary("/opt/perturbo.x") == "/opt/perturbo.x"
    assert runner.host_binary(Path("/ws"), "/opt/perturbo.x") is None


def test_runner_requires_dataset(monkeypatch, capsys):
    from fepdmc_toolkit.fepdmc import runner

    monkeypatch.delenv("FEPDMC_TOOLKIT_DATA", raising=False)
    assert runner.main(["ez", "lif_hole", "1"]) == 2
    assert "FEPDMC_TOOLKIT_DATA" in capsys.readouterr().err


def test_runner_reads_the_site_from_the_environment(monkeypatch):
    from fepdmc_toolkit.fepdmc import runner

    monkeypatch.setenv("FEPDMC_TOOLKIT_WS", "/w")
    monkeypatch.setenv("FEPDMC_TOOLKIT_DATA", "/d")
    monkeypatch.setenv("FEPDMC_TOOLKIT_IMAGE", "img:1")
    a = runner.parser().parse_args(["ez", "lif_elec"])
    assert (a.ws, a.data, a.image) == ("/w", "/d", "img:1")
    a = runner.parser().parse_args(["ez", "lif_elec", "--ws", "/x"])
    assert a.ws == "/x"


def test_runner_dry_run_prints_the_commands_and_writes_nothing(tmp_path, capsys):
    from fepdmc_toolkit.fepdmc import runner

    ws, data = tmp_path / "ws", tmp_path / "data"
    argv = ["ez", "lif_hole", "2", "--ws", str(ws), "--data", str(data), "--wqfix", "--dry-run"]
    assert runner.main(argv) == 0
    lines = capsys.readouterr().out.splitlines()
    assert len(lines) == 2 and all(line.startswith("docker run --platform linux/amd64 --rm ") for line in lines)
    assert "-e FUTUREB_WQFIX=1" in lines[0] and "-e FUTUREB_EXTFIX=0" in lines[0]
    assert "-e FUTUREB_SEED=20260925" in lines[0] and f"-e FUTUREB_SEED={20260925 + 7919}" in lines[1]
    assert f"-v {ws}:/work:rw" in lines[0] and f"-v {data}:/data:ro" in lines[0]
    assert runner.main(["prepare-tables", "sto", "--ws", str(ws), "--data", str(data), "--dry-run"]) == 0
    assert "tabulate" not in capsys.readouterr().err
    assert not ws.exists() and not data.exists()


def test_runner_names_the_missing_piece_before_starting_docker(tmp_path, monkeypatch, capsys):
    from fepdmc_toolkit.fepdmc import runner

    def no_docker(args):
        raise AssertionError("docker must not be started")

    monkeypatch.setattr(runner.shutil, "which", lambda name: "/usr/bin/docker")
    monkeypatch.setattr(runner.subprocess, "call", no_docker)
    ws, data = tmp_path / "ws", tmp_path / "data"
    argv = ["ez", "lif_hole", "1", "--ws", str(ws), "--data", str(data)]
    assert runner.main(argv) == 2
    assert "dataset directory not found" in capsys.readouterr().err
    data.mkdir()
    assert runner.main(argv) == 2
    assert "perturbo.x: not found" in capsys.readouterr().err
    exe = ws / runner.QE_DIR / runner.DEFAULT_BUILD / "pert-src" / "perturbo.x"
    exe.parent.mkdir(parents=True)
    exe.touch()
    assert runner.main(argv) == 2
    assert "run prepare-tables first" in capsys.readouterr().err
    monkeypatch.setattr(runner.shutil, "which", lambda name: None)
    assert runner.main(argv) == 2
    assert "docker not found" in capsys.readouterr().err


def test_runner_starts_one_container_per_chain(tmp_path, monkeypatch, capsys):
    from fepdmc_toolkit.fepdmc import runner

    started = []
    monkeypatch.setattr(runner.shutil, "which", lambda name: "/usr/bin/docker")
    monkeypatch.setattr(runner.subprocess, "call", lambda args: started.append(args) or 0)
    ws, data = tmp_path / "ws", tmp_path / "data"
    data.mkdir()
    (ws / "tables" / "lif_elec" / "gkq-20").mkdir(parents=True)
    exe = ws / runner.QE_DIR / runner.DEFAULT_BUILD / "pert-src" / "perturbo.x"
    exe.parent.mkdir(parents=True)
    exe.touch()
    argv = ["ez", "lif_elec", "3", "--ws", str(ws), "--data", str(data), "--parallel", "1", "--tag", "t"]
    assert runner.main(argv) == 0
    assert "exit codes [0, 0, 0]" in capsys.readouterr().out
    assert len(started) == 3 and all(args[:2] == ["docker", "run"] for args in started)
    chain = ws / "t" / "lif_elec" / "chain02"
    assert sorted(p.name for p in chain.iterdir()) == ["diagMC.in", "pert.in", "run_meta.json", "temper.in"]
    meta = json.loads((chain / "run_meta.json").read_text())
    assert meta["seed"] == 20260925 + 2 * 7919 and meta["switches"]["FUTUREB_RB"] == "1"
    assert " calc_mode = 'diagmc-EZ'" in (chain / "pert.in").read_text()


# ---------------------------------------------------------------- command line


def test_cli_dispatch_and_help(capsys):
    import fepdmc_toolkit
    from fepdmc_toolkit.fepdmc import cli

    assert cli.main(["--help"]) == 0
    assert cli.main([]) == 2
    capsys.readouterr()
    assert cli.main(["unknown"]) == 2
    assert "unknown command 'unknown'" in capsys.readouterr().err
    assert cli.main(["--version"]) == 0
    assert capsys.readouterr().out.strip() == f"fepdmc-toolkit {fepdmc_toolkit.__version__}"
    assert set(cli.COMMANDS) == {"prepare", "run", "validate", "compare", "analyze"}
    for command in cli.COMMANDS:
        with pytest.raises(SystemExit) as stop:
            cli.main([command, "--help"]) if command != "analyze" else cli.main([command, "analyze_q", "--help"])
        assert stop.value.code == 0
    assert capsys.readouterr().out.count("usage: fepdmc-toolkit ") == len(cli.COMMANDS)


def test_analysis_modules_import_without_scipy(monkeypatch):
    # scipy is optional (fepdmc-toolkit[analysis]); importing a module must not need it.
    import importlib

    from fepdmc_toolkit.fepdmc.analysis import ANALYSES

    for mod in list(sys.modules):
        if mod.startswith(("fepdmc_toolkit.fepdmc.analysis.", "scipy.")) or mod == "scipy":
            monkeypatch.delitem(sys.modules, mod)
    monkeypatch.setitem(sys.modules, "scipy", None)
    for name in ANALYSES:
        mod = importlib.import_module(f"fepdmc_toolkit.fepdmc.analysis.{name}")
        assert callable(getattr(mod, "main", None)), name
