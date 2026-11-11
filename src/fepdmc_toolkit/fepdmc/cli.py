"""fepdmc-toolkit: command line of the FEP-DMC toolkit.

  fepdmc-toolkit prepare  <tree> [--profile fixes|research] [--dry-run]
  fepdmc-toolkit run      {prepare-tables,ez} <material> [n_chains] [options]
  fepdmc-toolkit validate [--build TREE] [--kernels native,pure,mixed] [...]
  fepdmc-toolkit compare  <runs_A> <runs_C> <material> <E_bare> <tau_max> [...]
  fepdmc-toolkit analyze  <analysis> [arguments]      (analyze --list names them)
  fepdmc-toolkit --version

Run "fepdmc-toolkit <command> --help" for the options of each command.
"""

from __future__ import annotations

import sys
from importlib import import_module

COMMANDS = {
    "prepare": "fepdmc_toolkit.fepdmc.prepare",
    "run": "fepdmc_toolkit.fepdmc.runner",
    "validate": "fepdmc_toolkit.fepdmc.validate",
    "compare": "fepdmc_toolkit.fepdmc.pooled",
    "analyze": "fepdmc_toolkit.fepdmc.analysis",
}


def main(argv: list[str] | None = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    if not argv:
        print(__doc__)
        return 2
    command = argv[0]
    if command in ("-h", "--help"):
        print(__doc__)
        return 0
    if command in ("-V", "--version"):
        from fepdmc_toolkit import __version__

        print(f"fepdmc-toolkit {__version__}")
        return 0
    if command not in COMMANDS:
        print(f"fepdmc-toolkit: unknown command {command!r}\n", file=sys.stderr)
        print(__doc__, file=sys.stderr)
        return 2
    return int(import_module(COMMANDS[command]).main(argv[1:]) or 0)


if __name__ == "__main__":
    sys.exit(main())
