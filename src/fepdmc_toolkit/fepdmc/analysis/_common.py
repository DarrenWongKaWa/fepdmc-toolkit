"""Command-line plumbing shared by the analysis modules."""

from __future__ import annotations

import argparse
import json
import sys

# tau_max = 1/kT in 1/eV at T = 50 K, the temperature of the published runs.
TAU_MAX_50K = 232.09011


class AnalysisError(RuntimeError):
    """The run directories do not contain what the analysis needs."""


def make_parser(module: str, doc: str | None) -> argparse.ArgumentParser:
    """Argument parser of one analysis; ``module`` is its ``__name__``."""
    name = module.rsplit(".", 1)[-1]
    return argparse.ArgumentParser(
        prog=f"fepdmc-toolkit analyze {name}", description=doc, formatter_class=argparse.RawDescriptionHelpFormatter
    )


def emit(compute) -> int:
    """Print the result of ``compute()`` as JSON; report unusable input as an error (exit status 2)."""
    try:
        result = compute()
    except ImportError:
        print("error: this analysis needs scipy (pip install 'fepdmc-toolkit[analysis]')", file=sys.stderr)
        return 2
    except (AnalysisError, FileNotFoundError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    json.dump(result, sys.stdout, indent=2)
    sys.stdout.write("\n")
    return 0
